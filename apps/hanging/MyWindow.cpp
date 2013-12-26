/*
 * Copyright (c) 2011-2013, Georgia Tech Research Corporation
 * All rights reserved.
 *
 * Author(s): Karen Liu <karenliu@cc.gatech.edu>
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

#include "apps/hanging/MyWindow.h"

#include <cstdio>

#include "dart/common/Timer.h"
#include "dart/math/Helpers.h"
#include "dart/collision/CollisionDetector.h"
#include "dart/constraint/ConstraintDynamics.h"
#include "dart/constraint/PointConstraint.h"
#include "dart/dynamics/BodyNode.h"
#include "dart/dynamics/GenCoord.h"
#include "dart/dynamics/Skeleton.h"
#include "dart/simulation/World.h"
#include "dart/gui/GLFuncs.h"

MyWindow::MyWindow(dart::simulation::World* _world)
  : SimWindow()
{
  SimWindow::setWorld(_world);

  mForce = Eigen::Vector3d::Zero();
  mImpulseDuration = 0;

  Eigen::Vector3d gravity(0.0, -9.81, 0.0);
  mWorld->setGravity(gravity);

  //
  std::vector<int> genCoordIds0;
  genCoordIds0.push_back(4);

  // default standing pose
  std::vector<int> genCoordIds1;
  genCoordIds1.push_back(1);
  genCoordIds1.push_back(6);  // left hip
  genCoordIds1.push_back(14); // left knee
  genCoordIds1.push_back(17); // left ankle
  genCoordIds1.push_back(9);  // right hip
  genCoordIds1.push_back(15); // right knee
  genCoordIds1.push_back(19); // right ankle
  genCoordIds1.push_back(13); // lower back
  genCoordIds1.push_back(29); // left shoulder
  genCoordIds1.push_back(32); // right shoulder

  Eigen::VectorXd initConfig0(1);
  initConfig0 << -2.9;
  Eigen::VectorXd initConfig1(10);
  initConfig1 << -0.1, 0.2, -0.5, 0.3, 0.2, -0.5, 0.3, -0.1, 0.5, -0.5;

  mWorld->getSkeleton(0)->setConfig(genCoordIds0, initConfig0);
  mWorld->getSkeleton(1)->setConfig(genCoordIds1, initConfig1);

  // create controller
  mController = new Controller(mWorld->getSkeleton(1),
                               mWorld->getConstraintHandler(),
                               mWorld->getTimeStep());

  for (int i = 0; i < mWorld->getSkeleton(1)->getNumGenCoords(); i++)
    mController->setDesiredDof(i, mController->getSkeleton()->getGenCoord(i)->get_q());

  // initialize constraint on the hand
  dart::dynamics::BodyNode* bd
      = mWorld->getSkeleton(1)->getBodyNode("h_hand_left");
  dart::constraint::PointConstraint* point1
      = new dart::constraint::PointConstraint(
          bd, bd->getLocalCOM(), bd->getWorldCOM(), 1);
  mWorld->getConstraintHandler()->addConstraint(point1);
  bd = mWorld->getSkeleton(1)->getBodyNode("h_hand_right");
  dart::constraint::PointConstraint* point2
      = new dart::constraint::PointConstraint(
          bd, bd->getLocalCOM(), bd->getWorldCOM(), 1);
  mWorld->getConstraintHandler()->addConstraint(point2);
}

MyWindow::~MyWindow()
{
}

void MyWindow::timeStepping()
{
  mWorld->getSkeleton(1)->getBodyNode("h_pelvis")->addExtForce(mForce);

  mController->setConstrForces(
        mWorld->getConstraintHandler()->getTotalConstraintForce(1));
  mController->computeTorques(mWorld->getSkeleton(1)->get_q(),
                              mWorld->getSkeleton(1)->get_dq());
  mWorld->getSkeleton(1)->setInternalForceVector(mController->getTorques());

  mWorld->step();

  // for perturbation test
  mImpulseDuration--;
  if (mImpulseDuration <= 0)
  {
    mImpulseDuration = 0;
    mForce.setZero();
  }
}

void MyWindow::drawSkels()
{
  for (unsigned int i = 0; i < mWorld->getNumSkeletons(); i++)
    mWorld->getSkeleton(i)->draw(mRI);

  // draw handholds
  mRI->setPenColor(Eigen::Vector3d(0.2, 0.2, 0.2));
  mRI->pushMatrix();
  glTranslated(0.0, -0.06, -0.52);
  mRI->drawEllipsoid(Eigen::Vector3d(0.1, 0.1, 0.1));
  mRI->popMatrix();
  mRI->setPenColor(Eigen::Vector3d(0.2, 0.2, 0.2));
  mRI->pushMatrix();
  glTranslated(0.0, -0.06, 0.52);
  mRI->drawEllipsoid(Eigen::Vector3d(0.1, 0.1, 0.1));
  mRI->popMatrix();

  // draw arrow
  if (mImpulseDuration > 0)
  {
    dart::dynamics::BodyNode* bodyNode
        = mWorld->getSkeleton(1)->getBodyNode("h_pelvis");
    Eigen::Vector3d poa = bodyNode->getWorldTransform().translation();
    Eigen::Vector3d start = poa - mForce / 10.0;
    double len = mForce.norm() / 10.0;
    dart::gui::drawArrow3D(start, mForce, len, 0.05, 0.1);
  }
}

void MyWindow::keyboard(unsigned char key, int x, int y)
{
  switch(key)
  {
    case ' ':  // use space key to play or stop the motion
      mSimulating = !mSimulating;
      if (mSimulating)
      {
        mPlay = false;
        glutTimerFunc( mDisplayTimeout, refreshTimer, 0);
      }
      break;
    case 's':  // simulate one frame
      if (!mPlay)
      {
        mForce = Eigen::Vector3d::Zero();
        mWorld->step();
        bake();
        glutPostRedisplay();
      }
      break;
    case '1':
      mForce[0] = 20;
      mImpulseDuration = 10.0;
      std::cout << "push forward" << std::endl;
      break;
    case '2':
      mForce[0] = -10;
      mImpulseDuration = 10.0;
      std::cout << "push backward" << std::endl;
      break;
    case '3':
      mForce[2] = 50;
      mImpulseDuration = 10.0;
      std::cout << "push right" << std::endl;
      break;
    case '4':
      mForce[2] = -50;
      mImpulseDuration = 10.0;
      std::cout << "push left" << std::endl;
      break;
    case 'p':  // playBack
      mPlay = !mPlay;
      if (mPlay)
      {
        mSimulating = false;
        glutTimerFunc( mDisplayTimeout, refreshTimer, 0);
      }
      break;
    case '[':  // step backward
      if (!mSimulating)
      {
        mPlayFrame--;
        if(mPlayFrame < 0)
          mPlayFrame = 0;
        glutPostRedisplay();
      }
      break;
    case ']':  // step forwardward
      if (!mSimulating)
      {
        mPlayFrame++;
        if(mPlayFrame >= mBakedStates.size())
          mPlayFrame = 0;
        glutPostRedisplay();
      }
      break;
    case 'v':  // show or hide markers
      mShowMarkers = !mShowMarkers;
      break;

    default:
      Win3D::keyboard(key,x,y);
  }
  glutPostRedisplay();
}
