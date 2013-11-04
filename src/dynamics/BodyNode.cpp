/*
 * Copyright (c) 2011, Georgia Tech Research Corporation
 * All rights reserved.
 *
 * Author(s): Sehoon Ha <sehoon.ha@gmail.com>
 *            Jeongseok Lee <jslee02@gmail.com>
 * Date: 05/14/2013
 *
 * Geoorgia Tech Graphics Lab and Humanoid Robotics Lab
 *
 * Directed by Prof. C. Karen Liu and Prof. Mike Stilman
 * <karenliu@cc.gatech.edu> <mstilman@cc.gatech.edu>
 *
 * This file is provided under the following "BSD-style" License:
 *   Redistribution and use in source and binary forms, with or
 *   without modification, are permitted provided that the following
 *   conditions are met:
 *   * Redistributions of source code must retain the above copyright
 *     notice, this list of conditions and the following disclaimer.
 *   * Redistributions in binary form must reproduce the above
 *     copyright notice, this list of conditions and the following
 *     disclaimer in the documentation and/or other materials provided
 *     with the distribution.
 *   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND
 *   CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES,
 *   INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
 *   MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 *   DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR
 *   CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 *   SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 *   LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF
 *   USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED
 *   AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
 *   LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
 *   ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 *   POSSIBILITY OF SUCH DAMAGE.
 */

#include "dynamics/BodyNode.h"

#include <algorithm>
#include <iostream>
#include <queue>
#include <stack>

#include "common/Console.h"
#include "math/Helpers.h"
#include "renderer/RenderInterface.h"
#include "dynamics/Joint.h"
#include "dynamics/Shape.h"
#include "dynamics/Skeleton.h"
#include "dynamics/Marker.h"

namespace dart {
namespace dynamics {

int BodyNode::msBodyNodeCount = 0;

BodyNode::BodyNode(const std::string& _name)
    : mSkelIndex(-1),
      mName(_name),
      mIsCollidable(true),
      mIsColliding(false),
      mSkeleton(NULL),
      mParentJoint(NULL),
      mParentBodyNode(NULL),
      mChildBodyNodes(std::vector<BodyNode*>(0)),
      mGravityMode(true),
      mCenterOfMass(Eigen::Vector3d::Zero()),
      mMass(1.0),
      mIxx(1.0),
      mIyy(1.0),
      mIzz(1.0),
      mIxy(0.0),
      mIxz(0.0),
      mIyz(0.0),
      mI(Eigen::Matrix6d::Identity()),
      mW(Eigen::Isometry3d::Identity()),
      mV(Eigen::Vector6d::Zero()),
      mEta(Eigen::Vector6d::Zero()),
      mdV(Eigen::Vector6d::Zero()),
      mF(Eigen::Vector6d::Zero()),
      mFext(Eigen::Vector6d::Zero()),
      mFgravity(Eigen::Vector6d::Zero()),
      mAI(Eigen::Matrix6d::Identity()),
      mB(Eigen::Vector6d::Zero()),
      mBeta(Eigen::Vector6d::Zero()),
      mY(Eigen::Matrix6d::Zero()),
      mZ(Eigen::Matrix6d::Zero()),
      mID(BodyNode::msBodyNodeCount++)
{
}

BodyNode::~BodyNode()
{
    for (std::vector<Shape*>::const_iterator it = mVizShapes.begin();
         it != mVizShapes.end(); ++it)
        delete (*it);

    for (std::vector<Shape*>::const_iterator itColShape = mColShapes.begin();
         itColShape != mColShapes.end(); ++itColShape)
        if (mVizShapes.end() == std::find(mVizShapes.begin(), mVizShapes.end(), *itColShape))
            delete (*itColShape);

    for (std::vector<Marker*>::const_iterator it = mMarkers.begin();
         it != mMarkers.end(); ++it)
        delete (*it);

    if (mParentJoint)
        delete mParentJoint;
}

void BodyNode::setName(const std::string& _name)
{
    mName = _name;
}

const std::string& BodyNode::getName() const
{
    return mName;
}

void BodyNode::setGravityMode(bool _gravityMode)
{
    mGravityMode = _gravityMode;
}

bool BodyNode::getGravityMode() const
{
    return mGravityMode;
}

bool BodyNode::isCollidable() const
{
    return mIsCollidable;
}

void BodyNode::setCollidable(bool _isCollidable)
{
    mIsCollidable = _isCollidable;
}

void BodyNode::setMass(double _mass)
{
    assert(_mass >= 0.0 && "Negative mass is not allowable.");
    mMass = _mass;
    _updateGeralizedInertia();
}

double BodyNode::getMass() const
{
    return mMass;
}

BodyNode*BodyNode::getParentBodyNode() const
{
    return mParentBodyNode;
}

void BodyNode::addChildBodyNode(BodyNode* _body)
{
    assert(_body != NULL);

    mChildBodyNodes.push_back(_body);
    _body->mParentBodyNode = this;
}

BodyNode* BodyNode::getChildBodyNode(int _idx) const
{
    assert(0 <= _idx && _idx < mChildBodyNodes.size());

    return mChildBodyNodes[_idx];
}

int BodyNode::getNumChildBodyNodes() const
{
    return mChildBodyNodes.size();
}

void BodyNode::addMarker(Marker* _marker)
{
    mMarkers.push_back(_marker);
}

int BodyNode::getNumMarkers() const
{
    return mMarkers.size();
}

Marker* BodyNode::getMarker(int _idx) const
{
    return mMarkers[_idx];
}

bool BodyNode::dependsOn(int _genCoordIndex) const
{
    return std::binary_search(mDependentGenCoordIndices.begin(),
                              mDependentGenCoordIndices.end(),
                              _genCoordIndex);
}

int BodyNode::getNumDependentGenCoords() const
{
    return mDependentGenCoordIndices.size();
}

int BodyNode::getDependentGenCoord(int _arrayIndex) const
{
    assert(0 <= _arrayIndex && _arrayIndex < mDependentGenCoordIndices.size());
    return mDependentGenCoordIndices[_arrayIndex];
}

const Eigen::Isometry3d& BodyNode::getWorldTransform() const
{
    return mW;
}

const Eigen::Vector6d& BodyNode::getBodyVelocity() const
{
    return mV;
}

Eigen::Vector6d BodyNode::getWorldVelocity(const Eigen::Vector3d& _offset) const
{
    Eigen::Isometry3d T = mW;
    T.translation() = -_offset;
    return math::AdT(T, mV);
}

const Eigen::Vector6d&BodyNode::getBodyAcceleration() const
{
    return mdV;
}

Eigen::Vector6d BodyNode::getWorldAcceleration(
        const Eigen::Vector3d& _offset) const
{
    Eigen::Isometry3d T = mW;
    T.translation() = -_offset;
    return math::AdT(T, mdV);
}

const math::Jacobian&BodyNode::getBodyJacobian() const
{
    return mBodyJacobian;
}

math::Jacobian BodyNode::getWorldJacobian(const Eigen::Vector3d& _offset) const
{
    Eigen::Isometry3d T = mW;
    T.translation() = -_offset;
    return math::AdTJac(T, mBodyJacobian);
}

const math::Jacobian& BodyNode::getBodyJacobianTimeDeriv() const
{
    return mBodyJacobianTimeDeriv;
}

math::Jacobian BodyNode::getWorldJacobianTimeDeriv(
        const Eigen::Vector3d& _offset) const
{
    Eigen::Isometry3d T = mW;
    T.translation() = -_offset;
    return math::AdTJac(T, mBodyJacobianTimeDeriv);
}

void BodyNode::setColliding(bool _isColliding)
{
    mIsColliding = _isColliding;
}

bool BodyNode::isColliding()
{
    return mIsColliding;
}

void BodyNode::init(Skeleton* _skeleton, int _skeletonIndex)
{
    assert(_skeleton);

    mSkeleton = _skeleton;
    mSkelIndex = _skeletonIndex;
    mParentJoint->mSkelIndex = _skeletonIndex;

    //--------------------------------------------------------------------------
    // Fill the list of generalized coordinates this node depends on, and sort
    // it.
    //--------------------------------------------------------------------------
    if (mParentBodyNode)
        mDependentGenCoordIndices = mParentBodyNode->mDependentGenCoordIndices;
    else
        mDependentGenCoordIndices.clear();
    for (int i = 0; i < mParentJoint->getNumGenCoords(); i++)
        mDependentGenCoordIndices.push_back(mParentJoint->getGenCoord(i)->getSkeletonIndex());
    std::sort(mDependentGenCoordIndices.begin(), mDependentGenCoordIndices.end());

#ifndef NDEBUG
    // Check whether there is duplicated indices.
    for (int i = 0; i < mDependentGenCoordIndices.size() - 1; i++)
    {
        for (int j = i + 1; j < mDependentGenCoordIndices.size(); j++)
        {
            assert(mDependentGenCoordIndices[i] !=
                    mDependentGenCoordIndices[j] &&
                    "Duplicated index is found in mDependentGenCoordIndices.");
        }
    }
#endif

    //--------------------------------------------------------------------------
    // Fill the list of decendants of this body node.
    //--------------------------------------------------------------------------
    mDescendantBodyNodes.clear();

    // A. DFS (Depth First Search).
//    std::stack<BodyNode*> stack;
//    for (std::vector<BodyNode*>::const_reverse_iterator it =
//         mChildBodyNodes.rbegin(); it != mChildBodyNodes.rend(); ++it)
//    {
//        stack.push(*it);
//    }
//    while (!stack.empty())
//    {
//        BodyNode* itBodyNode = stack.top();
//        stack.pop();
//        mDescendantBodyNodes.push_back(itBodyNode);

//        for (int i = itBodyNode->getNumChildBodyNodes() - 1; i > -1 ; -i)
//        {
//            stack.push(itBodyNode->getChildBodyNode(i));
//        }
//    }

    // B. BFS (Breadth First Search)
    std::queue<BodyNode*> queue;
    for (std::vector<BodyNode*>::const_iterator it =
         mChildBodyNodes.begin(); it != mChildBodyNodes.end(); ++it)
    {
        queue.push(*it);
    }
    while (!queue.empty())
    {
        BodyNode* itBodyNode = queue.front();
        queue.pop();
        mDescendantBodyNodes.push_back(itBodyNode);

        for (int i = 0; i < itBodyNode->getNumChildBodyNodes(); ++i)
        {
            queue.push(itBodyNode->getChildBodyNode(i));
        }
    }

    //--------------------------------------------------------------------------
    // Fill all the nephew and younger brothers of this body
    //--------------------------------------------------------------------------
    mNephewAndYoungerBrotherBodyNodes.clear();
    // Add nephew of older brothers and younger brothers of this body to the
    // queue
    if (mParentBodyNode)
    {
        int i = 0;
        int numSibling = mParentBodyNode->getNumChildBodyNodes();
        for (; i < numSibling; ++i)
        {
            BodyNode* it = mParentBodyNode->getChildBodyNode(i);
            if (it == this)
            {
                ++i;
                break;
            }
            for (int j = 0; j < it->getNumChildBodyNodes(); ++j)
            {
                queue.push(it->getChildBodyNode(j));
            }
        }
        for (; i < numSibling; ++i)
        {
            BodyNode* it = mParentBodyNode->getChildBodyNode(i);
            queue.push(it);
        }
    }
    // Add all the nephew and younger brothers of this body
    while (!queue.empty())
    {
        BodyNode* it = queue.front();
        queue.pop();
        mNephewAndYoungerBrotherBodyNodes.push_back(it);

        for (int i = 0; i < it->getNumChildBodyNodes(); ++i)
        {
            queue.push(it->getChildBodyNode(i));
        }
    }

    for (std::vector<BodyNode*>::iterator it1 = mDescendantBodyNodes.begin();
         it1 != mDescendantBodyNodes.end(); ++it1)
    {
        assert(std::find(mNephewAndYoungerBrotherBodyNodes.begin(), mNephewAndYoungerBrotherBodyNodes.end(), *it1) == mNephewAndYoungerBrotherBodyNodes.end());
    }

    //--------------------------------------------------------------------------
    // Set dimensions of dynamics matrices and vectors.
    //--------------------------------------------------------------------------
    int numDepGenCoords = getNumDependentGenCoords();
    mBodyJacobian.setZero(6, numDepGenCoords);
    mBodyJacobianTimeDeriv.setZero(6, numDepGenCoords);
    mM.setZero(numDepGenCoords, numDepGenCoords);

    //--------------------------------------------------------------------------
    // Set dimensions of cache data for recursive algorithms
    //--------------------------------------------------------------------------
    int dof = mParentJoint->getNumGenCoords();
    mAI_S.setZero(6, dof);
    mPsi.setZero(dof, dof);
    mPsiK.setZero(dof, dof);
    mAlpha.setZero(dof);
    mP.setZero(dof, dof);
    mZ.setZero(6, dof);
    mA.setZero();
    mY.setIdentity();
}

void BodyNode::draw(renderer::RenderInterface* _ri,
                    const Eigen::Vector4d& _color,
                    bool _useDefaultColor,
                    int _depth) const
{
    if (_ri == NULL)
        return;

    _ri->pushMatrix();

    // render the self geometry
    mParentJoint->applyGLTransform(_ri);

    _ri->pushName((unsigned)mID);
    for(int i = 0; i < mVizShapes.size(); i++)
    {
        _ri->pushMatrix();
        mVizShapes[i]->draw(_ri, _color, _useDefaultColor);
        _ri->popMatrix();
    }
    _ri->popName();

    // render the subtree
    for (unsigned int i = 0; i < mChildBodyNodes.size(); i++)
    {
        mChildBodyNodes[i]->draw(_ri, _color, _useDefaultColor);
    }

    _ri->popMatrix();
}

void BodyNode::drawMarkers(renderer::RenderInterface* _ri,
                           const Eigen::Vector4d& _color,
                           bool _useDefaultColor) const
{
    if (!_ri)
        return;

    _ri->pushMatrix();

    mParentJoint->applyGLTransform(_ri);

    // render the corresponding mMarkerss
    for (unsigned int i = 0; i < mMarkers.size(); i++)
        mMarkers[i]->draw(_ri, true, _color, _useDefaultColor);

    for (unsigned int i = 0; i < mChildBodyNodes.size(); i++)
        mChildBodyNodes[i]->drawMarkers(_ri,_color, _useDefaultColor);

    _ri->popMatrix();
}

void BodyNode::updateTransform(bool _updateJacobian)
{
    mParentJoint->updateTransform();
    if (mParentBodyNode)
    {
        mW = mParentBodyNode->getWorldTransform()
             * mParentJoint->getLocalTransform();
    }
    else
    {
        mW = mParentJoint->getLocalTransform();
    }
    assert(math::verifyTransform(mW));

//    if (_updateJacobian == false)
//        return;

    mParentJoint->updateJacobian();

    //--------------------------------------------------------------------------
    // Jacobian update
    //
    // J = | J1 J2 ... Jn |
    //   = | Ad(T(i,i-1), J_parent) J_local |
    //
    //   J_parent: (6 x parentDOF)
    //    J_local: (6 x localDOF)
    //         Ji: (6 x 1) se3
    //          n: number of dependent coordinates
    //--------------------------------------------------------------------------

    const int numLocalDOFs = mParentJoint->getNumGenCoords();
    const int numParentDOFs = getNumDependentGenCoords()-numLocalDOFs;

    // Parent Jacobian
    if (mParentBodyNode != NULL)
    {
        assert(mParentBodyNode->mBodyJacobian.cols() + mParentJoint->getNumGenCoords()
               == mBodyJacobian.cols());

        for (int i = 0; i < numParentDOFs; ++i)
        {
            assert(mParentJoint);
            mBodyJacobian.col(i) = math::AdInvT(
                                       mParentJoint->getLocalTransform(),
                                       mParentBodyNode->mBodyJacobian.col(i));
        }
    }

    // Local Jacobian
    for(int i = 0; i < numLocalDOFs; i++)
    {
        mBodyJacobian.col(numParentDOFs + i) =
                mParentJoint->getLocalJacobian().col(i);
    }
}

void BodyNode::updateVelocity()
{
    //--------------------------------------------------------------------------
    // Body velocity update
    //
    // V(i) = Ad(T(i, i-1), V(i-1)) + S * dq
    //--------------------------------------------------------------------------

    if (mParentJoint->getNumGenCoords() > 0)
    {
        if (mParentBodyNode)
        {
            mV = math::AdInvT(mParentJoint->getLocalTransform(),
                              mParentBodyNode->getBodyVelocity()) +
                    mParentJoint->getLocalJacobian() * mParentJoint->get_dq();
        }
        else
        {
            mV = mParentJoint->getLocalJacobian() * mParentJoint->get_dq();
        }
    }

    assert(!math::isNan(mV));
}

void BodyNode::updateEta(bool _updateJacobianDeriv)
{
    mParentJoint->updateJacobianTimeDeriv();

    if (mParentJoint->getNumGenCoords() > 0)
    {
        mEta = math::ad(mV, mParentJoint->getLocalJacobian() *
                            mParentJoint->get_dq()) +
                            mParentJoint->getLocalJacobianTimeDeriv() *
                            mParentJoint->get_dq();

        assert(!math::isNan(mEta));
    }

    if (_updateJacobianDeriv == false)
        return;

    //--------------------------------------------------------------------------
    // Jacobian first derivative update
    //
    // dJ = | dJ1 dJ2 ... dJn |
    //   = | Ad(T(i,i-1), dJ_parent) dJ_local |
    //
    //   dJ_parent: (6 x parentDOF)
    //    dJ_local: (6 x localDOF)
    //         dJi: (6 x 1) se3
    //          n: number of dependent coordinates
    //--------------------------------------------------------------------------

    const int numLocalDOFs = mParentJoint->getNumGenCoords();
    const int numParentDOFs = getNumDependentGenCoords() - numLocalDOFs;

    // Parent Jacobian
    if (mParentBodyNode != NULL)
    {
        assert(mParentBodyNode->mBodyJacobianTimeDeriv.cols() + mParentJoint->getNumGenCoords()
               == mBodyJacobianTimeDeriv.cols());

        for (int i = 0; i < numParentDOFs; ++i)
        {
            assert(mParentJoint);
            Eigen::Vector6d dJi = math::AdInvT(mParentJoint->getLocalTransform(),
                                         mParentBodyNode->mBodyJacobianTimeDeriv.col(i));
            mBodyJacobianTimeDeriv.col(i) = dJi;
        }
    }

    // Local Jacobian
    for(int i = 0; i < numLocalDOFs; i++)
    {
        mBodyJacobianTimeDeriv.col(numParentDOFs + i) =
                mParentJoint->getLocalJacobianTimeDeriv().col(i);
    }
}

void BodyNode::updateAcceleration()
{
    // dV(i) = Ad(T(i, i-1), dV(i-1))
    //         + ad(V(i), S * dq) + dS * dq
    //         + S * ddq
    //       = Ad(T(i, i-1), dV(i-1))
    //         + eta
    //         + S * ddq

    if (mParentJoint->getNumGenCoords() > 0)
    {
        if (mParentBodyNode)
        {
            mdV = math::AdInvT(mParentJoint->getLocalTransform(),
                               mParentBodyNode->getBodyAcceleration()) +
                               mEta + mParentJoint->getLocalJacobian() *
                               mParentJoint->get_ddq();
        }
        else
        {
            mdV = mEta +
                  mParentJoint->getLocalJacobian() * mParentJoint->get_ddq();
        }
    }

    assert(!math::isNan(mdV));
}

void BodyNode::setInertia(double _Ixx, double _Iyy, double _Izz,
                          double _Ixy, double _Ixz, double _Iyz)
{
    assert(_Ixx >= 0.0);
    assert(_Iyy >= 0.0);
    assert(_Izz >= 0.0);

    mIxx = _Ixx;
    mIyy = _Iyy;
    mIzz = _Izz;

    mIxy = _Ixy;
    mIxz = _Ixz;
    mIyz = _Iyz;

    _updateGeralizedInertia();
}

void BodyNode::setLocalCOM(const Eigen::Vector3d& _com)
{
    mCenterOfMass = _com;

    _updateGeralizedInertia();
}

const Eigen::Vector3d& BodyNode::getLocalCOM() const
{
    return mCenterOfMass;
}

Eigen::Vector3d BodyNode::getWorldCOM() const
{
    return mW.linear() * mCenterOfMass;
}

Eigen::Matrix6d BodyNode::getInertia() const
{
    return mI;
}

int BodyNode::getSkeletonIndex() const
{
    return mSkelIndex;
}

void BodyNode::addVisualizationShape(Shape* _p)
{
    mVizShapes.push_back(_p);
}

int BodyNode::getNumVisualizationShapes() const
{
    return mVizShapes.size();
}

Shape*BodyNode::getVisualizationShape(int _idx) const
{
    return mVizShapes[_idx];
}

void BodyNode::addCollisionShape(Shape* _p)
{
    mColShapes.push_back(_p);
}

int BodyNode::getNumCollisionShapes() const
{
    return mColShapes.size();
}

Shape*BodyNode::getCollisionShape(int _idx) const
{
    return mColShapes[_idx];
}

Skeleton*BodyNode::getSkeleton() const
{
    return mSkeleton;
}

void BodyNode::setParentJoint(Joint* _joint)
{
    mParentJoint = _joint;
}

Joint*BodyNode::getParentJoint() const
{
    return mParentJoint;
}

void BodyNode::addExtForce(const Eigen::Vector3d& _offset,
                           const Eigen::Vector3d& _force,
                           bool _isOffsetLocal, bool _isForceLocal)
{
    Eigen::Isometry3d T = Eigen::Isometry3d::Identity();
    Eigen::Vector6d F = Eigen::Vector6d::Zero();

    if (_isOffsetLocal)
        T.translation() = _offset;
    else
        T.translation() = getWorldTransform().inverse() * _offset;

    if (_isForceLocal)
        F.tail<3>() = _force;
    else
        F.tail<3>() = mW.linear().transpose() * _force;

    mFext += math::dAdInvT(T, F);
}

void BodyNode::setExtForce(const Eigen::Vector3d& _offset,
                           const Eigen::Vector3d& _force,
                           bool _isOffsetLocal, bool _isForceLocal)
{
    Eigen::Isometry3d T = Eigen::Isometry3d::Identity();
    Eigen::Vector6d F = Eigen::Vector6d::Zero();

    if (_isOffsetLocal)
        T.translation() = _offset;
    else
        T.translation() = getWorldTransform().inverse() * _offset;

    if (_isForceLocal)
        F.tail<3>() = _force;
    else
        F.tail<3>() = mW.linear().transpose() * _force;

    mFext = math::dAdInvT(T, F);
}

void BodyNode::addExtTorque(const Eigen::Vector3d& _torque, bool _isLocal)
{
    if (_isLocal)
        mFext.head<3>() += _torque;
    else
        mFext.head<3>() += mW.linear() * _torque;
}

void BodyNode::setExtTorque(const Eigen::Vector3d& _torque, bool _isLocal)
{
    if (_isLocal)
        mFext.head<3>() = _torque;
    else
        mFext.head<3>() = mW.linear() * _torque;
}

const Eigen::Vector6d& BodyNode::getExternalForceLocal() const
{
    return mFext;
}

Eigen::Vector6d BodyNode::getExternalForceGlobal() const
{
    return math::dAdInvT(mW, mFext);
}

void BodyNode::addContactForce(const Eigen::Vector6d& _contactForce)
{
    mContactForces.push_back(_contactForce);
}

int BodyNode::getNumContactForces() const
{
    return mContactForces.size();
}

const Eigen::Vector6d& BodyNode::getContactForce(int _idx)
{
    assert(0 <= _idx && _idx < mContactForces.size());

    return mContactForces[_idx];
}

void BodyNode::clearContactForces()
{
    mContactForces.clear();
}

const Eigen::Vector6d& BodyNode::getBodyForce() const
{
    return mF;
}

double BodyNode::getKineticEnergy() const
{
    return 0.5 * mV.dot(mI * mV);
}

Eigen::Vector3d BodyNode::getLinearMomentum() const
{
    return (mI * mV).tail<3>();
}

Eigen::Vector3d BodyNode::getAngularMomentum(const Eigen::Vector3d& _pivot)
{
    Eigen::Isometry3d T = Eigen::Isometry3d::Identity();
    T.translation() = _pivot;
    return math::dAdT(T, mI * mV).head<3>();
}

void BodyNode::updateBodyForce(const Eigen::Vector3d& _gravity,
                               bool _withExternalForces)
{
    if (mGravityMode == true)
        mFgravity = mI * math::AdInvRLinear(mW, _gravity);
    else
        mFgravity.setZero();

    mF = mI * mdV;                // Inertial force
    if (_withExternalForces)
        mF -= mFext;              // External force
    mF -= mFgravity;              // Gravity force
    mF -= math::dad(mV, mI * mV); // Coriolis force

    for (std::vector<BodyNode*>::iterator iChildBody = mChildBodyNodes.begin();
         iChildBody != mChildBodyNodes.end();
         ++iChildBody)
    {
        Joint* childJoint = (*iChildBody)->getParentJoint();
        assert(childJoint != NULL);

        mF += math::dAdInvT(childJoint->getLocalTransform(),
                            (*iChildBody)->getBodyForce());
    }

    assert(!math::isNan(mF));
}

void BodyNode::updateGeneralizedForce(bool _withDampingForces)
{
    assert(mParentJoint != NULL);

    const math::Jacobian& J = mParentJoint->getLocalJacobian();

//    if (_withDampingForces)
//        mF -= mFDamp;

    assert(!math::isNan(J.transpose()*mF));

    mParentJoint->set_tau(J.transpose()*mF);
}

void BodyNode::updateArticulatedInertia(double _timeStep)
{
    assert(mParentJoint != NULL);

    // Articulated inertia
    mAI = mI;
    for (std::vector<BodyNode*>::const_iterator it = mChildBodyNodes.begin();
         it != mChildBodyNodes.end(); ++it)
    {
        mAI += math::transformInertia(
                    (*it)->getParentJoint()->getLocalTransform().inverse(),
                    (*it)->mPi);
    }
    assert(!math::isNan(mAI));

    // Cache data: PsiK and Psi
    mAI_S.noalias() = mAI * mParentJoint->getLocalJacobian();
    int dof = mParentJoint->getNumGenCoords();
    if (dof > 0)
    {
        Eigen::MatrixXd K = Eigen::MatrixXd::Zero(dof, dof);
        for (int i = 0; i < dof; ++i)
            K(i, i) = mParentJoint->getDampingCoefficient(i);
        mPsiK = (mParentJoint->getLocalJacobian().transpose() * mAI_S +
                 _timeStep * K).inverse();
        mPsi = (mParentJoint->getLocalJacobian().transpose() * mAI_S).inverse();
    }
    assert(!math::isNan(mPsiK));
    assert(!math::isNan(mPsi));

    // Cache data: Pi
    mPi = mAI;
    if (dof > 0)
        mPi.noalias() -= mAI_S*mPsiK*mAI_S.transpose();
    assert(!math::isNan(mPi));
}

void BodyNode::updateBiasForce(const Eigen::Vector3d& _gravity)
{
    // Bias force
    if (mGravityMode == true)
        mFgravity = mI * math::AdInvRLinear(mW, Eigen::Vector3d(_gravity));
    else
        mFgravity.setZero();
    mB = -math::dad(mV, mI*mV) - mFext - mFgravity;
    for (int i = 0; i < mContactForces.size(); ++i)
        mB -= mContactForces[i];
    for (std::vector<BodyNode*>::const_iterator it = mChildBodyNodes.begin();
         it != mChildBodyNodes.end(); ++it)
    {
        mB += math::dAdInvT((*it)->getParentJoint()->getLocalTransform(),
                            (*it)->mBeta);
    }
    assert(!math::isNan(mB));

    // Cache data: alpha
    int dof = mParentJoint->getNumGenCoords();
    if (dof > 0)
    {
        mAlpha = mParentJoint->get_tau() + mParentJoint->getDampingForces();
        for (int i = 0; i < dof; i++)
            mAlpha(i) += mSkeleton->getConstraintForceVector()[mParentJoint->getGenCoord(i)->getSkeletonIndex()];
        mAlpha -= mParentJoint->getLocalJacobian().transpose()*(mAI*mEta + mB);
        assert(!math::isNan(mAlpha));
    }

    // Cache data: beta
    mBeta = mB;
    if (dof > 0)
        mBeta += mAI*(mEta + mParentJoint->getLocalJacobian() * mPsiK * mAlpha);
    else
        mBeta += mAI*mEta;
    assert(!math::isNan(mBeta));
}

void BodyNode::update_ddq()
{
    if (mParentJoint->getNumGenCoords() == 0)
        return;

    Eigen::VectorXd ddq;
    if (mParentBodyNode)
    {
        ddq.noalias() = mPsiK*
                        (mAlpha -
                         mParentJoint->getLocalJacobian().transpose() * mAI *
                         math::AdInvT(mParentJoint->getLocalTransform(),
                                      mParentBodyNode->getBodyAcceleration())
                         );
    }
    else
    {
        ddq.noalias() = mPsiK*mAlpha;
    }

    mParentJoint->set_ddq(ddq);
    assert(!math::isNan(ddq));
}

void BodyNode::update_F_fs()
{
    mF.noalias() = mAI*mdV + mB;
    assert(!math::isNan(mF));
}

void BodyNode::updateZ()
{
    mZ.noalias() = math::dAdT(mParentJoint->getLocalTransform().inverse()) *
                   (mAI * mParentJoint->getLocalJacobian() * mPsi);
    assert(!math::isNan(mF));
}

void BodyNode::aggregateExternalForces(Eigen::VectorXd& _Fext)
{
    assert(mParentJoint != NULL);

    Eigen::VectorXd localForce = mBodyJacobian.transpose() * mFext;

    for(int i = 0; i < getNumDependentGenCoords(); i++)
        _Fext(mDependentGenCoordIndices[i]) += localForce(i);
}

void BodyNode::aggregateMassMatrix_OLD(Eigen::MatrixXd& _M)
{
    mM.triangularView<Eigen::Upper>() = mBodyJacobian.transpose() *
                                        mI *
                                        mBodyJacobian;
    mM.triangularView<Eigen::StrictlyLower>() = mM.transpose();

    for(int i = 0; i < getNumDependentGenCoords(); i++)
        for(int j = 0; j < getNumDependentGenCoords(); j++)
            _M(mDependentGenCoordIndices[i], mDependentGenCoordIndices[j]) += mM(i, j);
}

void BodyNode::aggregateMassMatrix(Eigen::MatrixXd& _M)
{
    BodyNode* currBodyNode = this;
    BodyNode* childBodyNode = NULL;

    int n = mParentJoint->getNumGenCoords();

    Eigen::MatrixXd M(n,n);

    const math::Jacobian& J = mParentJoint->getLocalJacobian();

    //--------------------------------------------------------------------------
    // M(j,j)
    //--------------------------------------------------------------------------
    math::Inertia tmpI = mI;

    for (int i = 0; i < getNumChildBodyNodes(); i++)
    {
        const Joint* childJoint         = mChildBodyNodes[i]->getParentJoint();
        const Eigen::Isometry3d& childT = childJoint->getLocalTransform();
        tmpI += mChildBodyNodes[i]->mIt * math::AdT(childT.inverse());
    }
    M = J.transpose() * tmpI * J;

    // Assigning
    for (int i = 0; i < n; i++)
    {
        int iSkelIdx = mParentJoint->getGenCoord(i)->getSkeletonIndex();
        for(int j = 0; j < n; j++)
        {
            int jSkelIdx = mParentJoint->getGenCoord(j)->getSkeletonIndex();
            _M(iSkelIdx, jSkelIdx) += M(i, j);
        }
    }

    childBodyNode = currBodyNode;
    currBodyNode = currBodyNode->getParentBodyNode();

    if (currBodyNode)
    {
        const Joint* childJoint         = childBodyNode->getParentJoint();
        const Eigen::Isometry3d& childT = childJoint->getLocalTransform();
        mIt = math::dAdT(childT.inverse()) * tmpI;
        tmpI = mIt;

        int m = currBodyNode->getParentJoint()->getNumGenCoords();
        Eigen::MatrixXd M2(m,n);
        const math::Jacobian& J2 = currBodyNode->getParentJoint()->getLocalJacobian();

        M2 = J2.transpose() * tmpI * J;

        // Assigning
        for(int i = 0; i < m; i++)
        {
            int iSkelIdx = currBodyNode->getParentJoint()->getGenCoord(i)->getSkeletonIndex();

            for(int j = 0; j < n; j++)
            {
                int jSkelIdx = mParentJoint->getGenCoord(j)->getSkeletonIndex();

                _M(iSkelIdx, jSkelIdx) += M2(i, j);
                _M(jSkelIdx, iSkelIdx) += M2(i, j);
            }
        }

        childBodyNode = currBodyNode;
        currBodyNode = currBodyNode->getParentBodyNode();
    }

    //--------------------------------------------------------------------------
    // M(j,k), k = lambda(j), lambda(lambda(j)), ... , base
    //--------------------------------------------------------------------------
    while (currBodyNode)
    {
        currBodyNode->mIt = math::Inertia::Zero();

        const Joint* childJoint         = childBodyNode->getParentJoint();
        const Eigen::Isometry3d& childT = childJoint->getLocalTransform();
        tmpI = math::dAdT(childT.inverse()) * tmpI;

        int m = currBodyNode->getParentJoint()->getNumGenCoords();
        Eigen::MatrixXd M2(m,n);
        const math::Jacobian& J2 = currBodyNode->getParentJoint()->getLocalJacobian();

        M2 = J2.transpose() * tmpI * J;

        // Assigning
        for(int i = 0; i < m; i++)
        {
            int iSkelIdx = currBodyNode->getParentJoint()->getGenCoord(i)->getSkeletonIndex();

            for(int j = 0; j < n; j++)
            {
                int jSkelIdx = mParentJoint->getGenCoord(j)->getSkeletonIndex();

                _M(iSkelIdx, jSkelIdx) += M2(i, j);
                _M(jSkelIdx, iSkelIdx) += M2(i, j);
            }
        }

        childBodyNode = currBodyNode;
        currBodyNode = currBodyNode->getParentBodyNode();
    }
}

void BodyNode::update_O_P_Z()
{
    mO = mAI * mParentJoint->getLocalJacobian() * mPsi;
    mP = math::Inertia::Identity() - mO * mParentJoint->getLocalJacobian().transpose();
    mZ = math::dAdT(mParentJoint->getLocalTransform().inverse()) * mO;
}

void BodyNode::aggregateInvMassMatrix(Eigen::MatrixXd& _MInv)
{
    int nGenCoords = mParentJoint->getNumGenCoords();
    const math::Jacobian& S = mParentJoint->getLocalJacobian();
    const Eigen::Isometry3d& T = mParentJoint->getLocalTransform();
    math::Inertia tmpInertia;

    // MInv(i, i) =
    //   Psi(i) + Z_{\lambda(i),i}^T A_{\lambda(i),\lambda(i)} Z_{\lambda(i),i}
    Eigen::MatrixXd MInv = mPsi;
    if (mParentBodyNode)
    {
        if (mParentBodyNode->getChildBodyNode(0) == this)
        {
            if (mParentBodyNode->mParentBodyNode)
            {
                //mParentBodyNode->mA = mParentBodyNode->mParentJoint  mParentBodyNode->mA;
            }
            else
            {
                mParentBodyNode->mA.setZero();
            }

            mParentBodyNode->mA +=
                    mParentBodyNode->mParentJoint->getLocalJacobian() *
                    mParentBodyNode->mPsi *
                    mParentBodyNode->mParentJoint->getLocalJacobian().transpose();
        }
        MInv += mZ.transpose() * mParentBodyNode->mA * mZ;
    }

    if (mChildBodyNodes.size() > 0)
    {
        // A_{\lambda(i),i} = A_{\lambda(i),\lambda(i)} Y_{\lambda(i),i}
        if (mParentBodyNode)
        {
            mA = mParentBodyNode->mA * mY;
            tmpInertia = mY.transpose() * mA;
        }
//        else
//        {
//            mA = math::Inertia::Zero();
//        }

        // Y
        mY = math::Inertia::Identity();
    }
    // Assigning
    for(int i = 0; i < nGenCoords; i++)
    {
        int iSkelIdx = mParentJoint->getGenCoord(i)->getSkeletonIndex();
        for(int j = 0; j < nGenCoords; j++)
        {
            int jSkelIdx = mParentJoint->getGenCoord(j)->getSkeletonIndex();
            _MInv(iSkelIdx, jSkelIdx) = MInv(i, j);
        }
    }

    // Descendants
    int nDescendants = mDescendantBodyNodes.size();
    for (int iDes = 0; iDes < nDescendants; ++iDes)
    {
        BodyNode* it = mDescendantBodyNodes[iDes];
        BodyNode* itParent = it->getParentBodyNode();
        Eigen::MatrixXd MInv2 = -mPsi * S.transpose() * itParent->mY * it->mZ;
        if (mParentBodyNode)
            MInv2 += mZ.transpose() * itParent->mA * it->mZ;
        if (it->getNumChildBodyNodes() > 0)
        {
            math::Inertia tmpY = math::dAdT(it->getParentJoint()->getLocalTransform().inverse()) * it->mP;

            // Y
            it->mY = itParent->mY * tmpY;

            // A
            if (it->mParentBodyNode)
                it->mA = itParent->mA * tmpY;
        }

        // Assigning
        int m = it->getParentJoint()->getNumGenCoords();
        for(int i = 0; i < nGenCoords; i++)
        {
            int iSkelIdx = mParentJoint->getGenCoord(i)->getSkeletonIndex();
            for(int j = 0; j < m; j++)
            {
                int jSkelIdx = it->getParentJoint()->getGenCoord(j)->getSkeletonIndex();
                _MInv(iSkelIdx, jSkelIdx) = MInv2(i, j);
                _MInv(jSkelIdx, iSkelIdx) = MInv2(i, j);
            }
        }
    }

    // Relatives
    int nRelatives = mNephewAndYoungerBrotherBodyNodes.size();
    for (int iRel = 0; iRel < nRelatives; ++iRel)
    {
        BodyNode* it = mNephewAndYoungerBrotherBodyNodes[iRel];
        BodyNode* itParent = it->getParentBodyNode();
        int m = it->getParentJoint()->getNumGenCoords();
        Eigen::MatrixXd MInv2 = Eigen::MatrixXd::Zero(nGenCoords,m);
        if (mParentBodyNode)
            MInv2 += mZ.transpose() * itParent->mA * it->mZ;
        if (it->getNumChildBodyNodes() > 0)
        {
            math::Inertia tmpY = math::dAdT(it->getParentJoint()->getLocalTransform().inverse()) * it->mP;

//            // Y
//            it->mY = itParent->mY * tmpY;

            // A
            if (it->mParentBodyNode)
                it->mA = itParent->mA * tmpY;
        }

        // Assigning
        for(int i = 0; i < nGenCoords; i++)
        {
            int iSkelIdx = mParentJoint->getGenCoord(i)->getSkeletonIndex();
            for(int j = 0; j < m; j++)
            {
                // TODO: need block asigning for performance
                int jSkelIdx = it->getParentJoint()->getGenCoord(j)->getSkeletonIndex();
                _MInv(iSkelIdx, jSkelIdx) = MInv2(i, j);
                _MInv(jSkelIdx, iSkelIdx) = MInv2(i, j);
            }
        }
    }

    // A
//    if (mChildBodyNodes.size() > 0)
//    {
//        mA = S * mPsi * S.transpose();
//        if (mParentBodyNode)
//            mA += tmpInertia;
//    }
}

void BodyNode::updateW()
{
    if (mParentJoint->getNumGenCoords() > 0)
    {
        if (mParentBodyNode)
        {
            mW2 = math::AdInvT(mParentJoint->getLocalTransform(),
                              mParentBodyNode->mW2) + mEta;
        }
        else
        {
            mW2 = mEta;
        }
    }
}

void BodyNode::aggregateCombinedVector(Eigen::VectorXd& _Cg, const Eigen::Vector3d& _gravity)
{
    //    mH =
    if (mGravityMode == true)
        mFgravity = mI * math::AdInvRLinear(mW, _gravity);
    else
        mFgravity.setZero();

    mH = mI * mdV;                // Inertial force
    //    if (_withExternalForces)
    //        mF -= mFext;              // External force
    mH -= mFgravity;              // Gravity force
    mH -= math::dad(mV, mI * mV); // Coriolis force

    for (std::vector<BodyNode*>::iterator iChildBody = mChildBodyNodes.begin();
         iChildBody != mChildBodyNodes.end();
         ++iChildBody)
    {
        Joint* childJoint = (*iChildBody)->getParentJoint();
        assert(childJoint != NULL);

        mH += math::dAdInvT(childJoint->getLocalTransform(), (*iChildBody)->mH);
    }
    Eigen::VectorXd Cg = mParentJoint->getLocalJacobian().transpose() * mH;

    // Assigning
    int nGenCoords = mParentJoint->getNumGenCoords();
    for(int i = 0; i < nGenCoords; i++)
    {
        int iSkelIdx = mParentJoint->getGenCoord(i)->getSkeletonIndex();
        {
            // TODO: need block asigning for performance
            _Cg(iSkelIdx) = Cg(i);
        }
    }
}

void BodyNode::_updateGeralizedInertia()
{
    // G = | I - m * [r] * [r]   m * [r] |
    //     |          -m * [r]     m * 1 |

    // m * r
    double mr0 = mMass * mCenterOfMass[0];
    double mr1 = mMass * mCenterOfMass[1];
    double mr2 = mMass * mCenterOfMass[2];

    // m * [r] * [r]
    double mr0r0 = mr0 * mCenterOfMass[0];
    double mr1r1 = mr1 * mCenterOfMass[1];
    double mr2r2 = mr2 * mCenterOfMass[2];
    double mr0r1 = mr0 * mCenterOfMass[1];
    double mr1r2 = mr1 * mCenterOfMass[2];
    double mr2r0 = mr2 * mCenterOfMass[0];

    mI(0,0) =  mIxx + mr1r1 + mr2r2;   mI(0,1) =  mIxy - mr0r1;           mI(0,2) =  mIxz - mr2r0;           assert(mI(0,3) == 0.0);   mI(0,4) = -mr2;           mI(0,5) =  mr1;
                                       mI(1,1) =  mIyy + mr2r2 + mr0r0;   mI(1,2) =  mIyz - mr1r2;           mI(1,3) =  mr2;           assert(mI(1,4) == 0.0);   mI(1,5) = -mr0;
                                                                          mI(2,2) =  mIzz + mr0r0 + mr1r1;   mI(2,3) = -mr1;           mI(2,4) =  mr0;           assert(mI(2,5) == 0.0);
                                                                                                             mI(3,3) =  mMass;         assert(mI(3,4) == 0.0);   assert(mI(3,5) == 0.0);
                                                                                                                                       mI(4,4) =  mMass;         assert(mI(4,5) == 0.0);
                                                                                                                                                                 mI(5,5) =  mMass;

    mI.triangularView<Eigen::StrictlyLower>() = mI.transpose();
}

void BodyNode::clearExternalForces()
{
    mFext.setZero();
    mContactForces.clear();
}

} // namespace dynamics
} // namespace dart

