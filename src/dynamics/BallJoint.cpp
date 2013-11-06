/*
 * Copyright (c) 2011, Georgia Tech Research Corporation
 * All rights reserved.
 *
 * Author(s): Jeongseok Lee <jslee02@gmail.com>
 * Date: 05/21/2013
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

#include "BallJoint.h"

#include "math/Helpers.h"
#include "math/Geometry.h"

namespace dart {
namespace dynamics {

BallJoint::BallJoint(const std::string& _name)
    : Joint(BALL, _name)
{
    mGenCoords.push_back(&mCoordinate[0]);
    mGenCoords.push_back(&mCoordinate[1]);
    mGenCoords.push_back(&mCoordinate[2]);

    mS = Eigen::Matrix<double,6,3>::Zero();
    mdS = Eigen::Matrix<double,6,3>::Zero();

    mDampingCoefficient.resize(3, 0);
}

BallJoint::~BallJoint()
{
}

inline void BallJoint::updateTransform()
{
    Eigen::Vector3d q(mCoordinate[0].get_q(),
                      mCoordinate[1].get_q(),
                      mCoordinate[2].get_q());

    mT = mT_ParentBodyToJoint *
            math::expAngular(q) *
            mT_ChildBodyToJoint.inverse();
}

inline void BallJoint::updateJacobian()
{
    Eigen::Vector3d q(mCoordinate[0].get_q(),
                      mCoordinate[1].get_q(),
                      mCoordinate[2].get_q());

    Eigen::Matrix3d J = math::expMapJac(q);

    Eigen::Vector6d J0;
    Eigen::Vector6d J1;
    Eigen::Vector6d J2;

    J0 << J(0,0), J(0,1), J(0,2), 0, 0, 0;
    J1 << J(1,0), J(1,1), J(1,2), 0, 0, 0;
    J2 << J(2,0), J(2,1), J(2,2), 0, 0, 0;

    mS.col(0) = math::AdT(mT_ChildBodyToJoint, J0);
    mS.col(1) = math::AdT(mT_ChildBodyToJoint, J1);
    mS.col(2) = math::AdT(mT_ChildBodyToJoint, J2);
}

inline void BallJoint::updateJacobianTimeDeriv()
{
    Eigen::Vector3d q(mCoordinate[0].get_q(),
                      mCoordinate[1].get_q(),
                      mCoordinate[2].get_q());
    Eigen::Vector3d dq(mCoordinate[0].get_dq(),
                       mCoordinate[1].get_dq(),
                       mCoordinate[2].get_dq());

    Eigen::Matrix3d dJ = math::expMapJacDot(q, dq);

    Eigen::Vector6d dJ0;
    Eigen::Vector6d dJ1;
    Eigen::Vector6d dJ2;

    dJ0 << dJ(0,0), dJ(0,1), dJ(0,2), 0, 0, 0;
    dJ1 << dJ(1,0), dJ(1,1), dJ(1,2), 0, 0, 0;
    dJ2 << dJ(2,0), dJ(2,1), dJ(2,2), 0, 0, 0;

    mdS.col(0) = math::AdT(mT_ChildBodyToJoint, dJ0);
    mdS.col(1) = math::AdT(mT_ChildBodyToJoint, dJ1);
    mdS.col(2) = math::AdT(mT_ChildBodyToJoint, dJ2);
}

void BallJoint::clampRotation()
{
//    Eigen::Vector3d exmap = get_q();

//    double theta = exmap.norm();
//    if (theta > M_PI)
//    {
//        exmap.normalize();
//        exmap *= theta-2*M_PI;

//        //BodyNode* node = jnt->getChildNode();
//        int firstIndex = 0;
//        if (jnt->getParentNode() != NULL)
//            firstIndex = jnt->getParentNode()->getNumDependentDofs();

//        // extract the local Jw
//        Eigen::Matrix3d oldJwBody;
//        for (int j = 0; j < 3; j++)
//        {
//            // XXX do not use node->mTq here because it's not computed if the recursive algorithm is used; instead the derivative matrix is (re)computed explicitly.
//            Eigen::Matrix3d omegaSkewSymmetric = jnt->getDeriv(jnt->getDof(j)).topLeftCorner<3,3>()
//                                          * node->getLocalTransform().topLeftCorner<3,3>().transpose();
//            oldJwBody.col(j) = dart_math::fromSkewSymmetric(omegaSkewSymmetric);
//        }

//        // the new Jw
//        Eigen::Matrix3d newJwBody;
//        // set new dof values to joint for derivative evaluation
//        for (int j = 0; j < 3; j++)
//            jnt->getDof(j)->setValue(exmap[j]);

//        // compute the new Jw from Rq*trans(R)
//        Eigen::Matrix4d Tbody = jnt->getLocalTransform();
//        for (int j = 0; j < 3; j++)
//        {
//            Eigen::Matrix3d omegaSkewSymmetric = jnt->getDeriv(jnt->getDof(j)).topLeftCorner<3,3>()*Tbody.topLeftCorner<3,3>().transpose();
//            newJwBody.col(j) = dart_math::fromSkewSymmetric(omegaSkewSymmetric);
//        }

//        // solve new_qdot s.t. newJw*new_qdot = w
//        // new_qdot = newJw.inverse()*w
//        Eigen::Vector3d old_qdot; // extract old_qdot
//        for (int j = 0; j < 3; j++)
//            old_qdot[j] = _qdot[jnt->getDof(j)->getSkelIndex()];
//        Eigen::Vector3d new_qdot = newJwBody.inverse()*oldJwBody*old_qdot;
//        for (int j = 0; j < 3; j++)
//        {
//            int dofIndex = jnt->getDof(j)->getSkelIndex();
//            _q[dofIndex] = exmap[j];
//            _qdot[dofIndex] = new_qdot[j];
//        }
//    }
}

} // namespace dynamics
} // namespace dart
