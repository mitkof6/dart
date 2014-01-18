/*
 * Copyright (c) 2014, Georgia Tech Research Corporation
 * All rights reserved.
 *
 * Author(s): Jeongseok Lee <jslee02@gmail.com>
 *
 * Georgia Tech Graphics Lab and Humanoid Robotics Lab
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
#include <iomanip>
#include <vector>

#include "dart/common/Console.h"
#include "dart/common/Timer.h"
#include "dart/math/Geometry.h"
#include "dart/math/Helpers.h"
#include "dart/dynamics/BodyNode.h"
#include "dart/dynamics/BoxShape.h"
#include "dart/dynamics/BallJoint.h"
#include "dart/dynamics/FreeJoint.h"
#include "dart/dynamics/PrismaticJoint.h"
#include "dart/dynamics/RevoluteJoint.h"
#include "dart/dynamics/Joint.h"
#include "dart/dynamics/TranslationalJoint.h"
#include "dart/dynamics/UniversalJoint.h"
#include "dart/dynamics/WeldJoint.h"
#include "dart/dynamics/EulerJoint.h"
#include "dart/dynamics/ScrewJoint.h"
#include "dart/dynamics/Skeleton.h"
#include "dart/simulation/World.h"
#include "dart/utils/Paths.h"
#include "dart/utils/SkelParser.h"

using namespace Eigen;

using namespace dart;
using namespace dynamics;

enum TypeOfDOF
{
  DOF_X, DOF_Y, DOF_Z, DOF_ROLL, DOF_PITCH, DOF_YAW
};

/// Add an end-effector to the last link of the given robot
void addEndEffector(Skeleton* robot, BodyNode* parent_node, Vector3d dim)
{
  // Create the end-effector node with a random dimension
  BodyNode* node = new BodyNode("ee");
  WeldJoint* joint = new WeldJoint("eeJoint");
  Eigen::Isometry3d T = Eigen::Isometry3d::Identity();
  T.translate(Eigen::Vector3d(0.0, 0.0, dim(2)));
  joint->setTransformFromParentBodyNode(T);
  Shape* shape = new BoxShape(Vector3d(0.2, 0.2, 0.2));
  node->setLocalCOM(Vector3d(0.0, 0.0, 0.0));
  node->setMass(1.0);
  node->addVisualizationShape(shape);
  node->addCollisionShape(shape);
  node->setParentJoint(joint);
  parent_node->addChildBodyNode(node);
  robot->addBodyNode(node);
}

/// Add a DOF to a given joint
Joint* create1DOFJoint(double val, double min, double max, int type)
{
  // Create the transformation based on the type
  Joint* newJoint = NULL;
  if(type == DOF_X)
    newJoint = new PrismaticJoint(Eigen::Vector3d(1.0, 0.0, 0.0));
  else if(type == DOF_Y)
    newJoint = new PrismaticJoint(Eigen::Vector3d(0.0, 1.0, 0.0));
  else if(type == DOF_Z)
    newJoint = new PrismaticJoint(Eigen::Vector3d(0.0, 0.0, 1.0));
  else if(type == DOF_YAW)
    newJoint = new RevoluteJoint(Eigen::Vector3d(0.0, 0.0, 1.0));
  else if(type == DOF_PITCH)
    newJoint = new RevoluteJoint(Eigen::Vector3d(0.0, 1.0, 0.0));
  else if(type == DOF_ROLL)
    newJoint = new RevoluteJoint(Eigen::Vector3d(1.0, 0.0, 0.0));
  // Add the transformation to the joint, set the min/max values and set it to
  // the skeleton
  newJoint->getGenCoord(0)->set_q(val);
  newJoint->getGenCoord(0)->set_qMin(min);
  newJoint->getGenCoord(0)->set_qMax(max);

  return newJoint;
}

/// Creates a N link manipulator with the given dimensions where the first link
/// rotates around z-axis and second rotates around x in the zero configuration.
Skeleton* createNLinkRobot(int _n, Vector3d dim, TypeOfDOF type,
                           bool finished = false)
{
  assert(_n > 0);

  Skeleton* robot = new Skeleton();
  robot->setSelfCollidable(false);
  double mass = 1.0;

  // Create the first link, the joint with the ground and its shape
  BodyNode* parent_node = NULL;
  BodyNode* node = new BodyNode("link1");
  Joint* joint = create1DOFJoint(0.0, -DART_PI, DART_PI, type);
  joint->setName("joint1");
  joint->setDampingCoefficient(0, 0.01);
  Shape* shape = new BoxShape(dim);
  node->setLocalCOM(Vector3d(0.0, 0.0, dim(2)/2.0));
  node->addVisualizationShape(shape);
  node->addCollisionShape(shape);
  node->setMass(mass);
  node->setParentJoint(joint);
  robot->addBodyNode(node);
  parent_node = node;

  // Create links iteratively
  for (int i = 1; i < _n; ++i)
  {
    std::ostringstream ssLink;
    std::ostringstream ssJoint;
    ssLink << "link" << i;
    ssLink << "joint" << i;
    node = new BodyNode(ssLink.str());
    joint = create1DOFJoint(0.0, -DART_PI, DART_PI, type);
    joint->setName(ssJoint.str());
    Eigen::Isometry3d T = Eigen::Isometry3d::Identity();
    T.translate(Eigen::Vector3d(0.0, 0.0, dim(2)));
    joint->setTransformFromParentBodyNode(T);
    joint->setDampingCoefficient(0, 0.01);
    shape = new BoxShape(dim);
    node->setLocalCOM(Vector3d(0.0, 0.0, dim(2)/2.0));
    node->addVisualizationShape(shape);
    node->addCollisionShape(shape);
    node->setMass(mass);
    node->setParentJoint(joint);
    parent_node->addChildBodyNode(node);
    robot->addBodyNode(node);
    parent_node = node;
  }

  // If finished, initialize the skeleton
  if(finished)
  {
    addEndEffector(robot, node, dim);
    robot->init();
  }
  return robot;
}

int main(int argc, char* argv[])
{
  common::Timer timer;
  std::vector<int> dofs;

  // Set 1
#ifdef NDEBUG // Release mode
  int numItr = 1000;
#else          // Debug mode
  int numItr = 10;
#endif
  dofs.push_back(5);
  dofs.push_back(10);
  dofs.push_back(15);
  dofs.push_back(20);
  dofs.push_back(25);
  dofs.push_back(30);
  dofs.push_back(35);
  dofs.push_back(40);
  dofs.push_back(45);
  dofs.push_back(50);

  // Set 2
  //#ifdef NDEBUG // Release mode
  //    int numItr = 500;
  //#else          // Debug mode
  //    int numItr = 5;
  //#endif
  //    dofs.push_back(10);
  //    dofs.push_back(20);
  //    dofs.push_back(30);
  //    dofs.push_back(40);
  //    dofs.push_back(50);
  //    dofs.push_back(60);
  //    dofs.push_back(70);
  //    dofs.push_back(80);
  //    dofs.push_back(90);
  //    dofs.push_back(100);

  // Set 3
  //#ifdef NDEBUG // Release mode
  //    int numItr = 5;
  //#else          // Debug mode
  //    int numItr = 1;
  //#endif
  //    dofs.push_back(100);
  //    dofs.push_back(200);
  //    dofs.push_back(300);
  //    dofs.push_back(400);
  //    dofs.push_back(500);
  //    dofs.push_back(600);
  //    dofs.push_back(700);

  std::vector<double> newResult(dofs.size(), 0);
  std::vector<double> newResultWithEOM(dofs.size(), 0);

  for (int i = 0; i < dofs.size(); ++i)
  {
    simulation::World world;
    dynamics::Skeleton* skeleton =
        createNLinkRobot(dofs[i], Eigen::Vector3d::Ones(), DOF_X, true);
    world.addSkeleton(skeleton);

    // Random state
    Eigen::VectorXd state = skeleton->getState();
    for (int k = 0; k < state.size(); ++k)
    {
      // TODO: The range is [-0.4pi, 0.4pi] until we resolve
      //       singular Jacobian issue.
      state[k] = math::random(-DART_PI*0.4, DART_PI*0.4);
    }
    skeleton->setState(state);

    timer.start();
    for (int j = 0; j < numItr; j++)
    {
      world.step();
    }
    timer.stop();
    newResult[i] = timer.getLastElapsedTime();
  }

  for (int i = 0; i < dofs.size(); ++i)
  {
    simulation::World world;
    dynamics::Skeleton* skeleton =
        createNLinkRobot(dofs[i], Eigen::Vector3d::Ones(), DOF_X, true);
    world.addSkeleton(skeleton);

    // Random state
    Eigen::VectorXd state = skeleton->getState();
    for (int k = 0; k < state.size(); ++k)
    {
      // TODO: The range is [-0.4pi, 0.4pi] until we resolve
      //       singular Jacobian issue.
      state[k] = math::random(-DART_PI*0.4, DART_PI*0.4);
    }
    skeleton->setState(state);

    timer.start();
    for (int j = 0; j < numItr; j++)
    {
      world.step();
      Eigen::MatrixXd M    = skeleton->getMassMatrix();
      Eigen::MatrixXd MInv = skeleton->getInvMassMatrix();
      Eigen::VectorXd Cg   = skeleton->getCombinedVector();
      Eigen::VectorXd C    = skeleton->getCoriolisForceVector();
      Eigen::VectorXd g    = skeleton->getGravityForceVector();
      Eigen::VectorXd Fext = skeleton->getExternalForceVector();
      Eigen::MatrixXd J    = skeleton->getBodyNode(skeleton->getNumBodyNodes()
                                                   - 1)->getBodyJacobian();
    }
    timer.stop();
    newResultWithEOM[i] = timer.getLastElapsedTime();
  }

  std::cout << "--------------------------------------------------------------"
            << std::endl;
  std::cout << std::setw(12) << " dof";
  for (int i = 0; i < dofs.size(); ++i)
    std::cout << " | " << std::setw(4) << dofs[i] << "";
  std::cout << std::endl;
  std::cout << "--------------------------------------------------------------"
            << std::endl;
  std::cout << std::setw(12) << " new w/o EOM";
  for (int i = 0; i < dofs.size(); ++i)
    std::cout << " | " << std::setw(4) << newResult[i] << "";
  std::cout << std::endl;
  std::cout << std::setw(12) << " new w EOM";
  for (int i = 0; i < dofs.size(); ++i)
    std::cout << " | " << std::setw(4) << newResultWithEOM[i] << "";
  std::cout << std::endl;

  return 0;
}
