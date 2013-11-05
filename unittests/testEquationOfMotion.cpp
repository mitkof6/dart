/*
 * Copyright (c) 2011, Georgia Tech Research Corporation
 * All rights reserved.
 *
 * Author(s): Jeongseok Lee <jslee02@gmail.com>
 * Date: 05/23/2013
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

#include <iostream>
#include <gtest/gtest.h>
#include "TestHelpers.h"

#include "math/Geometry.h"
#include "math/Helpers.h"
#include "dynamics/BallJoint.h"
#include "dynamics/FreeJoint.h"
#include "dynamics/PrismaticJoint.h"
#include "dynamics/RevoluteJoint.h"
#include "dynamics/Skeleton.h"
#include "dynamics/TranslationalJoint.h"
#include "dynamics/UniversalJoint.h"
#include "dynamics/WeldJoint.h"
#include "dynamics/EulerJoint.h"
#include "dynamics/ScrewJoint.h"
#include "simulation/World.h"
#include "utils/Paths.h"
#include "utils/SkelParser.h"

using namespace dart;
using namespace math;
using namespace dynamics;

#define EOM_TOL 0.01

/******************************************************************************/
class EOM : public testing::Test
{
public:
    void equationsOfMotionTest(const std::string& _fileName);
};

/******************************************************************************/
void EOM::equationsOfMotionTest(const std::string& _fileName)
{
    simulation::World* myWorld = utils::SkelParser::readSkelFile(_fileName);
    EXPECT_TRUE(myWorld != NULL);

    myWorld->setGravity(Eigen::Vector3d(0.0, 0.0, -9.81));
    myWorld->setTimeStep(0.001);


    double simTime = 0.001;
    double timeStep = myWorld->getTimeStep();
    int nSteps = simTime / timeStep;

    for (int i = 0; i < nSteps; i++)
    {
        myWorld->step();

        for (int j = 0; j < myWorld->getNumSkeletons(); ++j)
        {
            dynamics::Skeleton* skeleton = myWorld->getSkeleton(j);

            // mass matrix
            Eigen::MatrixXd M_OLD = skeleton->getMassMatrix_OLD();
            Eigen::MatrixXd M = skeleton->getMassMatrix();
            if (!equals(M_OLD, M))
            {
                std::cout << "M_OLD: " << M_OLD << std::endl;
                std::cout << "M    : " << M << std::endl;
            }
            EXPECT_TRUE(equals(M_OLD, M));

            // inverse mass matrix
            Eigen::MatrixXd MInv_OLD = skeleton->getInvMassMatrix_OLD();
            Eigen::MatrixXd MInv     = skeleton->getInvMassMatrix();
            if (!equals(MInv_OLD, MInv))
            {
                std::cout << "MInv_OLD: \n" << MInv_OLD << std::endl;
                std::cout << "MInv    : \n" << MInv << std::endl;
                std::cout << "Diff    : \n" << MInv - MInv_OLD << std::endl;
            }
            EXPECT_TRUE(equals(MInv_OLD, MInv));

            // Coriolis and gravity force
            Eigen::VectorXd Cg_OLD = skeleton->getCombinedVector_OLD();
            Eigen::VectorXd Cg     = skeleton->getCombinedVector();
            if (!equals(Cg_OLD, Cg))
            {
                std::cout << "Cg_OLD: \n" << Cg_OLD << std::endl;
                std::cout << "Cg    : \n" << Cg << std::endl;
            }
            EXPECT_TRUE(equals(Cg_OLD, Cg));
        }
    }

    int nRandomTest = 100;
    for (int i = 0; i < nRandomTest; i++)
    {
        for (int j = 0; j < myWorld->getNumSkeletons(); ++j)
        {
            dynamics::Skeleton* skeleton = myWorld->getSkeleton(j);

            Eigen::VectorXd state = skeleton->getState();
            for (int j = 0; j < state.size(); ++j)
                state[j] = math::random(-DART_PI, DART_PI);
            skeleton->setState(state);

            // mass matrix
            Eigen::MatrixXd M_OLD = skeleton->getMassMatrix_OLD();
            Eigen::MatrixXd M = skeleton->getMassMatrix();
            if (!equals(M_OLD, M))
            {
                std::cout << "M_OLD: " << M_OLD << std::endl;
                std::cout << "M    : " << M << std::endl;
            }
            EXPECT_TRUE(equals(M_OLD, M));

            // inverse mass matrix
//            Eigen::MatrixXd MInv_OLD = skeleton->getInvMassMatrix_OLD();
//            Eigen::MatrixXd MInv     = skeleton->getInvMassMatrix();
//            if (!equals(MInv_OLD, MInv))
//            {
//                std::cout << "MInv_OLD: \n" << MInv_OLD << std::endl;
//                std::cout << "MInv    : \n" << MInv << std::endl;
//            }
            //EXPECT_TRUE(equals(MInv_OLD, MInv));

            // Coriolis and gravity force
            Eigen::VectorXd Cg_OLD = skeleton->getCombinedVector_OLD();
            Eigen::VectorXd Cg     = skeleton->getCombinedVector();
            if (!equals(Cg_OLD, Cg))
            {
                std::cout << "Cg_OLD: \n" << Cg_OLD << std::endl;
                std::cout << "Cg    : \n" << Cg << std::endl;
            }
            EXPECT_TRUE(equals(Cg_OLD, Cg));
        }
    }
}

/******************************************************************************/
//TEST_F(EOM, SINGLE_PENDULUM)
//{
//    equationsOfMotionTest(DART_DATA_PATH"/skel/test/single_pendulum.skel");
//}

/******************************************************************************/
//TEST_F(EOM, DOUBLE_PENDULUM)
//{
//    equationsOfMotionTest(DART_DATA_PATH"/skel/test/double_pendulum.skel");
//}

/******************************************************************************/
//TEST_F(EOM, SIMPLE_TREE_STRUCTURE)
//{
//    equationsOfMotionTest(
//                DART_DATA_PATH"/skel/test/simple_tree_structure.skel");
//}

/******************************************************************************/
//TEST_F(EOM, TREE_STRUCTURE)
//{
//    equationsOfMotionTest(
//                DART_DATA_PATH"/skel/test/tree_structure.skel");
//}

/******************************************************************************/
TEST_F(EOM, FULL_BODY)
{
    equationsOfMotionTest(DART_DATA_PATH"/skel/fullbody1.skel");
}

/******************************************************************************/
int main(int argc, char* argv[])
{
	::testing::InitGoogleTest(&argc, argv);
	return RUN_ALL_TESTS();
}


