/*
 * Copyright (c) 2011-2014, Georgia Tech Research Corporation
 * All rights reserved.
 *
 * Author(s): Sumit Jain <sumit@cc.gatech.edu>,
 *            Jeongseok Lee <jslee02@gmail.com>
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

#include <Eigen/Dense>
#include <gtest/gtest.h>

//#include "TestHelpers.h"

//#include "dart/common/Console.h"
//#include "dart/math/Geometry.h"
//#include "dart/math/Helpers.h"
//#include "dart/dynamics/BodyNode.h"
//#include "dart/dynamics/Skeleton.h"
//#include "dart/simulation/World.h"
//#include "dart/utils/SkelParser.h"
//#include "dart/utils/Paths.h"
#include "dart/dynamics/State.h"

#include <boost/any.hpp>

//==============================================================================
TEST(StateTest, Any)
{
  using namespace Eigen;
  using namespace dart;
  using namespace dynamics;

//  Any n;
//  EXPECT_TRUE(n.isNull());

//  Isometry3d T1 = Isometry3d::Identity();
//  boost::any a1 = T1;
////  std::type_info ti = a1.type();
//  std::string name = a1.type().name();

//  std::cout << "name: " << name.c_str() << std::endl;
////  Isometry3d T1 = Isometry3d::Identity();
//  Any mya1;// = T1;
//  mya1 = T1;
////  EXPECT_TRUE(a1.not_null());
////  EXPECT_TRUE(a1.is<Isometry3d>());
////  EXPECT_TRUE(!a1.is<int>());
////  a = T;
}

//==============================================================================
int main(int argc, char* argv[]) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}

