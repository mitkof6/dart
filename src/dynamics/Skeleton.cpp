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

#include <queue>
#include <boost/math/special_functions/fpclassify.hpp>

#include "math/Geometry.h"
#include "math/Helpers.h"
#include "dynamics/BodyNode.h"
#include "dynamics/GenCoord.h"
#include "dynamics/Joint.h"
#include "dynamics/Marker.h"
#include "dynamics/Skeleton.h"

namespace dart {
namespace dynamics {

Skeleton::Skeleton(const std::string& _name)
    : GenCoordSystem(),
      mName(_name),
      mIsSelfCollidable(false),
      mTimeStep(0.001),
      mGravity(Eigen::Vector3d(0.0, 0.0, -9.81)),
      mTotalMass(0.0),
      mIsMobile(true),
      mIsMassMatrixDirty_OLD(true),
      mIsMassMatrixDirty(true),
      mIsMassInvMatrixDirty_OLD(true),
      mIsMassInvMatrixDirty(true),
      mIsCoriolisVectorDirty_OLD(true),
      mIsCoriolisVectorDirty(true),
      mIsGravityForceVectorDirty_OLD(true),
      mIsGravityForceVectorDirty(true),
      mIsCombinedVectorDirty_OLD(true),
      mIsCombinedVectorDirty(true),
      mIsExternalForceVectorDirty_OLD(true),
      mIsExternalForceVectorDirty(true),
      mIsDampingForceVectorDirty(true)
{
}

Skeleton::~Skeleton()
{
    for (std::vector<BodyNode*>::const_iterator it = mBodyNodes.begin();
         it != mBodyNodes.end(); ++it)
        delete (*it);
}

void Skeleton::setName(const std::string& _name)
{
    mName = _name;
}

const std::string& Skeleton::getName() const
{
    return mName;
}

void Skeleton::setSelfCollidable(bool _isSelfCollidable)
{
    mIsSelfCollidable = _isSelfCollidable;
}

bool Skeleton::isSelfCollidable() const
{
    return mIsSelfCollidable;
}

void Skeleton::setMobile(bool _isMobile)
{
    mIsMobile = _isMobile;
}

bool Skeleton::isMobile() const
{
    return mIsMobile;
}

void Skeleton::setTimeStep(double _timeStep)
{
    assert(_timeStep > 0.0);
    mTimeStep = _timeStep;
}

double Skeleton::getTimeStep() const
{
    return mTimeStep;
}

void Skeleton::setGravity(const Eigen::Vector3d& _gravity)
{
    mGravity = _gravity;
}

const Eigen::Vector3d& Skeleton::getGravity() const
{
    return mGravity;
}

double Skeleton::getMass() const
{
    return mTotalMass;
}

void Skeleton::init(double _timeStep, const Eigen::Vector3d& _gravity)
{
    // Set timestep and gravity
    setTimeStep(_timeStep);
    setGravity(_gravity);

    // Rearrange the list of body nodes with BFS (Breadth First Search)
    std::queue<BodyNode*> queue;
    queue.push(mBodyNodes[0]);
    mBodyNodes.clear();
    while (!queue.empty())
    {
        BodyNode* itBodyNode = queue.front();
        queue.pop();
        mBodyNodes.push_back(itBodyNode);
        for (int i = 0; i < itBodyNode->getNumChildBodyNodes(); ++i)
            queue.push(itBodyNode->getChildBodyNode(i));
    }

    // Initialize body nodes and generalized coordinates
    mGenCoords.clear();
    for(int i = 0; i < getNumBodyNodes(); ++i)
    {
        Joint* joint = mBodyNodes[i]->getParentJoint();
        for (int j = 0; j < joint->getNumGenCoords(); ++j)
        {
            joint->getGenCoord(j)->setSkeletonIndex(mGenCoords.size());
            mGenCoords.push_back(joint->getGenCoord(j));
        }
        mBodyNodes[i]->init(this, i);
        mBodyNodes[i]->updateTransform();
        mBodyNodes[i]->updateVelocity();
        mBodyNodes[i]->updateEta();
    }
    for (std::vector<BodyNode*>::reverse_iterator it = mBodyNodes.rbegin();
         it != mBodyNodes.rend(); ++it)
    {
        (*it)->updateArticulatedInertia(mTimeStep);
    }

//    for (int i = 0; i < mBodyNodes.size(); ++i)
//    {
//        std::cout << "BodyNode " << i << ":" << mBodyNodes[i]->mSkelIndex << ", " << mBodyNodes[i]->mName << std::endl;
//    }

    // Set dimension of dynamics quantities
    int dof = getNumGenCoords();
    mM_OLD = Eigen::MatrixXd::Zero(dof, dof);
    mM    = Eigen::MatrixXd::Zero(dof, dof);
    mMInv_OLD = Eigen::MatrixXd::Zero(dof, dof);
    mMInv = Eigen::MatrixXd::Zero(dof, dof);
    mCvec_OLD = Eigen::VectorXd::Zero(dof);
    mCvec = Eigen::VectorXd::Zero(dof);
    mG_OLD = Eigen::VectorXd::Zero(dof);
    mG    = Eigen::VectorXd::Zero(dof);
    mCg_OLD = Eigen::VectorXd::Zero(dof);
    mCg   = Eigen::VectorXd::Zero(dof);
    set_tau(Eigen::VectorXd::Zero(dof));
    mFext_OLD = Eigen::VectorXd::Zero(dof);
    mFext = Eigen::VectorXd::Zero(dof);
    mFc   = Eigen::VectorXd::Zero(dof);
    mFd   = Eigen::VectorXd::Zero(dof);

    // Calculate mass
    mTotalMass = 0.0;
    for(int i = 0; i < getNumBodyNodes(); i++)
        mTotalMass += getBodyNode(i)->getMass();
}

void Skeleton::addBodyNode(BodyNode* _body)
{
    assert(_body && _body->getParentJoint());

    mBodyNodes.push_back(_body);
}

int Skeleton::getNumBodyNodes() const
{
    return mBodyNodes.size();
}

BodyNode* Skeleton::getRootBodyNode() const
{
    // We assume that the first element of body nodes is root.
    return mBodyNodes[0];
}

BodyNode* Skeleton::getBodyNode(int _idx) const
{
    return mBodyNodes[_idx];
}

BodyNode* Skeleton::getBodyNode(const std::string& _name) const
{
    assert(!_name.empty());

    for (std::vector<BodyNode*>::const_iterator itrBody = mBodyNodes.begin();
         itrBody != mBodyNodes.end();
         ++itrBody) {
        if ((*itrBody)->getName() == _name)
            return *itrBody;
    }

    return NULL;
}

Joint* Skeleton::getJoint(int _idx) const
{
    return mBodyNodes[_idx]->getParentJoint();
}

Joint* Skeleton::getJoint(const std::string& _name) const
{
    assert(!_name.empty());

    for (std::vector<BodyNode*>::const_iterator it = mBodyNodes.begin();
         it != mBodyNodes.end();
         ++it)
    {
        if ((*it)->getParentJoint()->getName() == _name)
            return (*it)->getParentJoint();
    }

    return NULL;
}

Marker* Skeleton::getMarker(const std::string& _name) const
{
    assert(!_name.empty());

    for (std::vector<BodyNode*>::const_iterator it = mBodyNodes.begin();
         it != mBodyNodes.end(); ++it)
    {
        for (int i = 0; i < (*it)->getNumMarkers(); ++i)
        {
            if ((*it)->getMarker(i)->getName() == _name)
                return (*it)->getMarker(i);
        }
    }

    return NULL;
}

Eigen::VectorXd Skeleton::getConfig(const std::vector<int>& _id) const
{
    Eigen::VectorXd q(_id.size());

    for(unsigned int i = 0; i < _id.size(); i++)
        q[i] = mGenCoords[_id[i]]->get_q();

    return q;
}

Eigen::VectorXd Skeleton::getConfig() const
{
    return get_q();
}

void Skeleton::setConfig(const std::vector<int>& _id, const Eigen::VectorXd& _config)
{
    for( unsigned int i = 0; i < _id.size(); i++ )
        mGenCoords[_id[i]]->set_q(_config(i));

    for (std::vector<BodyNode*>::iterator itrBody = mBodyNodes.begin();
         itrBody != mBodyNodes.end(); ++itrBody)
    {
        (*itrBody)->updateTransform();
    }
}

void Skeleton::setConfig(const Eigen::VectorXd& _config)
{
    set_q(_config);

    for (std::vector<BodyNode*>::iterator itrBody = mBodyNodes.begin();
         itrBody != mBodyNodes.end(); ++itrBody)
    {
        (*itrBody)->updateTransform();
    }

    for (std::vector<BodyNode*>::reverse_iterator it = mBodyNodes.rbegin();
         it != mBodyNodes.rend(); ++it)
    {
        (*it)->updateArticulatedInertia(mTimeStep);
    }

    mIsMassMatrixDirty = true;
    mIsMassInvMatrixDirty = true;
    mIsCoriolisVectorDirty = true;
    mIsGravityForceVectorDirty = true;
    mIsCombinedVectorDirty = true;
    mIsExternalForceVectorDirty = true;
    //mIsDampingForceVectorDirty = true;

    mIsMassMatrixDirty_OLD = true;
    mIsMassInvMatrixDirty_OLD = true;
    mIsCoriolisVectorDirty_OLD = true;
    mIsGravityForceVectorDirty_OLD = true;
    mIsCombinedVectorDirty_OLD = true;
    mIsExternalForceVectorDirty_OLD = true;
//    mIsDampingForceVectorDirty_OLD = true;
}

void Skeleton::setState(const Eigen::VectorXd& _state)
{
    set_q(_state.head(_state.size() / 2));
    set_dq(_state.tail(_state.size() / 2));
    
    for (std::vector<BodyNode*>::iterator it = mBodyNodes.begin();
         it != mBodyNodes.end(); ++it)
    {
        (*it)->getParentJoint()->clampRotation();
        (*it)->updateTransform();
        (*it)->updateVelocity();
        (*it)->updateEta();
    }

    for (std::vector<BodyNode*>::reverse_iterator it = mBodyNodes.rbegin();
         it != mBodyNodes.rend(); ++it)
    {
        (*it)->updateArticulatedInertia(mTimeStep);
    }

    mIsMassMatrixDirty = true;
    mIsMassInvMatrixDirty = true;
    mIsCoriolisVectorDirty = true;
    mIsGravityForceVectorDirty = true;
    mIsCombinedVectorDirty = true;
    mIsExternalForceVectorDirty = true;
    mIsDampingForceVectorDirty = true;

    mIsMassMatrixDirty_OLD = true;
    mIsMassInvMatrixDirty_OLD = true;
    mIsCoriolisVectorDirty_OLD = true;
    mIsGravityForceVectorDirty_OLD = true;
    mIsCombinedVectorDirty_OLD = true;
    mIsExternalForceVectorDirty_OLD = true;
    //mIsDampingForceVectorDirty_OLD = true;
}

Eigen::VectorXd Skeleton::getState()
{
    Eigen::VectorXd state(2 * mGenCoords.size());
    state << get_q(), get_dq();
    return state;
}

/// Returns true if the two matrices are equal within the given bound
template <class MATRIX>
bool equals (const Eigen::DenseBase<MATRIX>& A, const Eigen::DenseBase<MATRIX>& B, double tol = 1e-5) {

    // Get the matrix sizes and sanity check the call
    const size_t n1 = A.cols(), m1 = A.rows();
    const size_t n2 = B.cols(), m2 = B.rows();
    if(m1!=m2 || n1!=n2) return false;

    // Check each index
    for(size_t i=0; i<m1; i++) {
        for(size_t j=0; j<n1; j++) {
            if(boost::math::isnan(A(i,j)) ^ boost::math::isnan(B(i,j)))
                return false;
            else if(fabs(A(i,j) - B(i,j)) > tol)
                return false;
        }
    }

    // If no problems, the two matrices are equal
    return true;
}

const Eigen::MatrixXd& Skeleton::getMassMatrix_OLD()
{
    if (mIsMassMatrixDirty_OLD)
        _updateMassMatrix_OLD();

    return mM_OLD;
}

const Eigen::MatrixXd& Skeleton::getMassMatrix()
{
    if (mIsMassMatrixDirty)
        _updateMassMatrix();

    Eigen::MatrixXd M = getMassMatrix_OLD();

    assert(equals(mM, M));

    return mM;
}

const Eigen::MatrixXd& Skeleton::getInvMassMatrix_OLD()
{
    if (mIsMassInvMatrixDirty_OLD)
        _updateInvMassMatrix_OLD();

    return mMInv_OLD;
}

const Eigen::MatrixXd& Skeleton::getInvMassMatrix()
{
    if (mIsMassInvMatrixDirty)
        _updateInvMassMatrix();

//    Eigen::MatrixXd MInv = getInvMassMatrix_OLD();

//    assert(equals(mMInv, MInv));

    return mMInv;
}

const Eigen::VectorXd& Skeleton::getCoriolisForceVector_OLD()
{
    if (mIsCoriolisVectorDirty_OLD)
        _updateCoriolisForceVector_OLD();

    return mCvec_OLD;
}

const Eigen::VectorXd& Skeleton::getCoriolisForceVector()
{
    if (mIsCoriolisVectorDirty)
        _updateCoriolisForceVector();

    return mCvec;
}

const Eigen::VectorXd& Skeleton::getGravityForceVector_OLD()
{
    if (mIsGravityForceVectorDirty_OLD)
        _updateGravityForceVector_OLD();

    return mG_OLD;
}

const Eigen::VectorXd& Skeleton::getGravityForceVector()
{
    if (mIsGravityForceVectorDirty)
        _updateGravityForceVector();

    return mG;
}

const Eigen::VectorXd& Skeleton::getCombinedVector_OLD()
{
    if (mIsCombinedVectorDirty_OLD)
        _updateCombinedVector_OLD();

    return mCg_OLD;
}

const Eigen::VectorXd& Skeleton::getCombinedVector()
{
    if (mIsCombinedVectorDirty)
        _updateCombinedVector();

    return mCg;
}

const Eigen::VectorXd& Skeleton::getExternalForceVector_OLD()
{
    if (mIsExternalForceVectorDirty_OLD)
        _updateExternalForceVector_OLD();

    return mFext_OLD;
}

const Eigen::VectorXd& Skeleton::getExternalForceVector()
{
    if (mIsExternalForceVectorDirty)
        _updateExternalForceVector();

    return mFext;
}

Eigen::VectorXd Skeleton::getInternalForceVector() const
{
    return get_tau();
}

const Eigen::VectorXd& Skeleton::getDampingForceVector()
{
    if (mIsDampingForceVectorDirty)
        _updateDampingForceVector();

    return mFd;
}

const Eigen::VectorXd& Skeleton::getConstraintForceVector()
{
    return mFc;
}

void Skeleton::draw(renderer::RenderInterface* _ri,
                    const Eigen::Vector4d& _color,
                    bool _useDefaultColor) const
{
    getRootBodyNode()->draw(_ri, _color, _useDefaultColor);
}

void Skeleton::drawMarkers(renderer::RenderInterface* _ri,
                           const Eigen::Vector4d& _color,
                           bool _useDefaultColor) const
{
    getRootBodyNode()->drawMarkers(_ri, _color, _useDefaultColor);
}

void Skeleton::_updateMassMatrix_OLD()
{
    assert(mM_OLD.cols() == getNumGenCoords() && mM_OLD.rows() == getNumGenCoords());
    assert(getNumGenCoords() > 0);

    // Calcualtion mass matrix, M
    mM_OLD.setZero();
    for (std::vector<BodyNode*>::iterator it = mBodyNodes.begin();
         it != mBodyNodes.end(); ++it)
    {
        (*it)->aggregateMassMatrix_OLD(mM_OLD);
    }

    mIsMassMatrixDirty_OLD = false;
}

void Skeleton::_updateMassMatrix()
{
    assert(mM.cols() == getNumGenCoords() && mM.rows() == getNumGenCoords());
    assert(getNumGenCoords() > 0);

    // Calcualtion mass matrix, M
    mM.setZero();
    for (std::vector<BodyNode*>::reverse_iterator it = mBodyNodes.rbegin();
         it != mBodyNodes.rend(); ++it)
    {
        (*it)->aggregateMassMatrix(mM);
        //(*it)->aggregateMassMatrix_OLD(mM);
    }

    mIsMassMatrixDirty = false;
}

void Skeleton::_updateInvMassMatrix_OLD()
{
    assert(mMInv_OLD.cols() == getNumGenCoords() &&
           mMInv_OLD.rows() == getNumGenCoords());
    assert(getNumGenCoords() > 0);

    // Update mass matrix
    if (mIsMassMatrixDirty_OLD)
        _updateMassMatrix_OLD();

    // Inverse of mass matrix
    int n = getNumGenCoords();
    mMInv_OLD = mM_OLD.ldlt().solve(Eigen::MatrixXd::Identity(n,n));

    mIsMassInvMatrixDirty_OLD = false;
}

void Skeleton::_updateInvMassMatrix()
{
    assert(mMInv.cols() == getNumGenCoords() &&
           mMInv.rows() == getNumGenCoords());
    assert(getNumGenCoords() > 0);

    // Prepare cache data
    for (std::vector<BodyNode*>::iterator it = mBodyNodes.begin();
         it != mBodyNodes.end(); ++it)
    {
        (*it)->update_O_P_Z();
    }

    // Inverse of mass matrix
    mMInv.setZero();
    for (std::vector<BodyNode*>::iterator it = mBodyNodes.begin();
         it != mBodyNodes.end(); ++it)
    {
        (*it)->aggregateInvMassMatrix(mMInv);
    }

    mIsMassInvMatrixDirty = false;
}

void Skeleton::_updateCoriolisForceVector_OLD()
{
    assert(mCvec_OLD.size() == getNumGenCoords());
    assert(getNumGenCoords() > 0);

    // Save current tau and gravity
    Eigen::VectorXd tau_old = get_tau();
    Eigen::Vector3d oldGravity = getGravity();

    // Set ddq as zero and gravity as zero
    set_ddq(Eigen::VectorXd::Zero(getNumGenCoords()));
    setGravity(Eigen::Vector3d::Zero());

    // M(q) * ddq + b(q,dq) = tau
    computeInverseDynamicsLinear(true);
    mCvec_OLD = get_tau();

    // Restore the torque
    set_tau(tau_old);
    setGravity(oldGravity);

    mIsCoriolisVectorDirty_OLD = false;
}

void Skeleton::_updateCoriolisForceVector()
{
    assert(mCvec.size() == getNumGenCoords());
    assert(getNumGenCoords() > 0);

    mCvec.setZero();
    for (std::vector<BodyNode*>::reverse_iterator it = mBodyNodes.rbegin();
         it != mBodyNodes.rend(); ++it)
    {
        (*it)->aggregateCoriolisForceVector(mCvec);
    }

    mIsCoriolisVectorDirty = false;
}

void Skeleton::_updateGravityForceVector_OLD()
{
    assert(mG_OLD.size() == getNumGenCoords());
    assert(getNumGenCoords() > 0);

    // Calcualtion mass matrix, M
    mG_OLD.setZero();
    for (std::vector<BodyNode*>::reverse_iterator it = mBodyNodes.rbegin();
         it != mBodyNodes.rend(); ++it)
    {
        (*it)->aggregateGravityForceVector_OLD(mG_OLD, mGravity);
    }

    mIsGravityForceVectorDirty_OLD = false;
}

void Skeleton::_updateGravityForceVector()
{
    assert(mG.size() == getNumGenCoords());
    assert(getNumGenCoords() > 0);

    // Calcualtion mass matrix, M
    mG.setZero();
    for (std::vector<BodyNode*>::reverse_iterator it = mBodyNodes.rbegin();
         it != mBodyNodes.rend(); ++it)
    {
        (*it)->aggregateGravityForceVector(mG, mGravity);
    }

    mIsGravityForceVectorDirty = false;
}

void Skeleton::_updateCombinedVector_OLD()
{
    assert(mCg_OLD.size() == getNumGenCoords());
    assert(getNumGenCoords() > 0);

    // Save current tau
    Eigen::VectorXd tau_old = get_tau();

    // Set ddq as zero
    set_ddq(Eigen::VectorXd::Zero(getNumGenCoords()));

    // M(q) * ddq + b(q,dq) = tau
    computeInverseDynamicsLinear(true);
    mCg_OLD = get_tau();

    // Restore the torque
    set_tau(tau_old);

    mIsCombinedVectorDirty_OLD = false;
}

void Skeleton::_updateCombinedVector()
{
    assert(mCg.size() == getNumGenCoords());
    assert(getNumGenCoords() > 0);

    mCg.setZero();
    for (std::vector<BodyNode*>::iterator it = mBodyNodes.begin();
         it != mBodyNodes.end(); ++it)
    {
        (*it)->updateW();
    }
    for (std::vector<BodyNode*>::reverse_iterator it = mBodyNodes.rbegin();
         it != mBodyNodes.rend(); ++it)
    {
        (*it)->aggregateCombinedVector(mCg, mGravity);
    }

    mIsCombinedVectorDirty = false;
}

void Skeleton::_updateExternalForceVector_OLD()
{
    assert(mFext_OLD.size() == getNumGenCoords());
    assert(getNumGenCoords() > 0);

    // Clear external force.
    mFext_OLD.setZero();
    for (std::vector<BodyNode*>::iterator itr = mBodyNodes.begin();
         itr != mBodyNodes.end(); ++itr)
        (*itr)->aggregateExternalForces_OLD(mFext_OLD);

    mIsExternalForceVectorDirty_OLD = false;
}

void Skeleton::_updateExternalForceVector()
{
    assert(mFext.size() == getNumGenCoords());
    assert(getNumGenCoords() > 0);

    // Clear external force.
    mFext.setZero();
    for (std::vector<BodyNode*>::reverse_iterator itr = mBodyNodes.rbegin();
         itr != mBodyNodes.rend(); ++itr)
        (*itr)->aggregateExternalForces(mFext);

    mIsExternalForceVectorDirty = false;
}

void Skeleton::_updateDampingForceVector()
{
    assert(mFd.size() == getNumGenCoords());
    assert(getNumGenCoords() > 0);

    // Clear external force.
    mFd.setZero();

    for (std::vector<BodyNode*>::iterator itr = mBodyNodes.begin();
         itr != mBodyNodes.end(); ++itr)
    {
        Eigen::VectorXd jointDampingForce = (*itr)->getParentJoint()->getDampingForces();
        for (int i = 0; i < jointDampingForce.size(); i++)
        {
            mFd((*itr)->getParentJoint()->getGenCoord(i)->getSkeletonIndex()) =
                    jointDampingForce(i);
        }
    }
}

void Skeleton::computeInverseDynamicsLinear(bool _computeJacobian,
                                            bool _computeJacobianDeriv,
                                            bool _withExternalForces,
                                            bool _withDampingForces)
{
    // Skip immobile or 0-dof skeleton
    if (getNumGenCoords() == 0)
        return;

    // Forward recursion
    for (std::vector<BodyNode*>::iterator it
         = mBodyNodes.begin(); it != mBodyNodes.end(); ++it)
    {
        (*it)->updateAcceleration();
    }

    // Backward recursion
    for (std::vector<BodyNode*>::reverse_iterator it
         = mBodyNodes.rbegin(); it != mBodyNodes.rend(); ++it)
    {
        (*it)->updateBodyForce(mGravity, _withExternalForces);
        (*it)->updateGeneralizedForce(_withDampingForces);
    }
}

void Skeleton::clearExternalForceVector()
{
    for (std::vector<BodyNode*>::iterator it = mBodyNodes.begin();
         it != mBodyNodes.end(); ++it)
    {
        (*it)->clearExternalForces();
    }
}

void Skeleton::computeForwardDynamicsID()
{
    // Skip immobile or 0-dof skeleton
    if (!isMobile() || getNumGenCoords() == 0)
        return;

    Eigen::VectorXd qddot = getInvMassMatrix()
                            * (-getCombinedVector()
                               + getExternalForceVector()
                               + getInternalForceVector()
                               + getDampingForceVector()
                               + getConstraintForceVector());

    this->set_ddq(qddot);
}

void Skeleton::computeForwardDynamicsFS()
{
    // Skip immobile or 0-dof skeleton
    if (!isMobile() || getNumGenCoords() == 0)
        return;

    // Backward recursion
    for (std::vector<BodyNode*>::reverse_iterator it = mBodyNodes.rbegin();
         it != mBodyNodes.rend(); ++it)
    {
        (*it)->updateBiasForce(mGravity);
    }

    // Forward recursion
    for (std::vector<BodyNode*>::iterator it = mBodyNodes.begin();
         it != mBodyNodes.end(); ++it)
    {
        (*it)->update_ddq();
        (*it)->updateAcceleration();
        (*it)->update_F_fs();
    }
}

void Skeleton::setInternalForceVector(const Eigen::VectorXd& _forces)
{
    set_tau(_forces);
}

void Skeleton::setMinInternalForceVector(const Eigen::VectorXd& _minForces)
{
    set_tauMin(_minForces);
}

Eigen::VectorXd Skeleton::getMinInternalForces() const
{
    return get_tauMin();
}

void Skeleton::setMaxInternalForceVector(const Eigen::VectorXd& _maxForces)
{
    set_tauMax(_maxForces);
}

Eigen::VectorXd Skeleton::getMaxInternalForceVector() const
{
    return get_tauMax();
}

void Skeleton::clearInternalForceVector()
{
    set_tau(Eigen::VectorXd::Zero(getNumGenCoords()));
}

void Skeleton::setConstraintForceVector(const Eigen::VectorXd& _Fc)
{
    mFc = _Fc;
}

Eigen::Vector3d Skeleton::getWorldCOM()
{
    Eigen::Vector3d com(0, 0, 0);

    assert(mTotalMass != 0);
    const int nNodes = getNumBodyNodes();

    for(int i = 0; i < nNodes; i++)
    {
        BodyNode* bodyNode = getBodyNode(i);
        com += bodyNode->getMass() *
               (bodyNode->getWorldTransform() * bodyNode->getLocalCOM());
    }

    return com / mTotalMass;
}

} // namespace dynamics
} // namespace dart
