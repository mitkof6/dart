/*
 * Copyright (c) 2011-2013, Georgia Tech Research Corporation
 * All rights reserved.
 *
 * Author(s): Jie (Jay) Tan <jtan34@cc.gatech.edu>
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

#include "dart/lcpsolver/LCPSolver.h"

#include <cstdio>

#include "dart/lcpsolver/Lemke.h"
#include "dart/lcpsolver/lcp.h"
#include "dart/lcpsolver/misc.h"

namespace dart {
namespace lcpsolver {

LCPSolver::LCPSolver()
{
}

LCPSolver::~LCPSolver()
{
}

bool LCPSolver::Solve(const Eigen::MatrixXd& _A,
                      const Eigen::VectorXd& _b,
                      Eigen::VectorXd* _x,
                      int _numContacts,
                      double _mu,
                      int _numDir,
                      bool _bUseODESolver)
{
  if (!_bUseODESolver)
  {
    int err = Lemke(_A, _b, _x);
    return (err == 0);
  }
  else
  {
    assert(_numDir >= 4);

    int n = _A.rows();    // A is square matrix
    int nSkip = dPAD(n);

    // Allocate intermediate variable for ODE
    double* A      = new double[n * nSkip];
    double* b      = new double[n];
    double* x      = new double[n];
    double* w      = new double[n]();  // initialize as 0
    double* lo     = new double[n]();  // initialize as 0
    double* hi     = new double[n];
    int*    findex = new int[n];

    // Fill values of A, b, w, lo, hi, and findex for normal contact forces
    for (int i = 0; i < n; ++i)
    {
      for (int j = 0; j < n; ++j)
        A[i * nSkip + j] = _A(i, j);
      b[i] = -_b[i];
      hi[i] = dInfinity;
      findex[i] = -1;
    }

    // Fill values of findex, lo, hi for directional contact forces
    for (int i = 0; i < _numContacts; ++i)
    {
      int idx0 = _numContacts + i * 2;
      int idx1 = idx0 + 1;

      findex[idx0] = i;
      findex[idx1] = i;

      lo[idx0] = -_mu;
      lo[idx1] = -_mu;

      hi[idx0] = _mu;
      hi[idx1] = _mu;
    }

    // dClearUpperTriangle(A,n);

    // Solve LCP
    dSolveLCP(n, A, x, b, w, 0, lo, hi, findex);

//    for (int i = 0; i < n; i++) {
//      if (w[i] < 0.0 && abs(x[i] - hi[i]) > 0.000001)
//        cout << "w[" << i << "] is negative, but x is " << x[i] << endl;
//      else if (w[i] > 0.0 && abs(x[i] - lo[i]) > 0.000001)
//        cout << "w[" << i << "] is positive, but x is " << x[i]
//                << " lo is " <<  lo[i] << endl;
//      else if (abs(w[i]) < 0.000001 && (x[i] > hi[i] || x[i] < lo[i]))
//        cout << "w[i] " << i << " is zero, but x is " << x[i] << endl;
//    }

    // Resize _x to n
    _x->setZero(n);

    // Assign the result of lcp into _x
    for (int i = 0; i < n; ++i)
      (*_x)[i] = x[i];

    // checkIfSolution(reducedA, reducedb, _x);

    // Release intermediate variables for ODE
    delete[] A;
    delete[] b;
    delete[] x;
    delete[] w;
    delete[] lo;
    delete[] hi;
    delete[] findex;

    return 1;
  }
}

void LCPSolver::transferToODEFormulation(const Eigen::MatrixXd& _A,
                                         const Eigen::VectorXd& _b,
                                         Eigen::MatrixXd* _AOut,
                                         Eigen::VectorXd* _bOut,
                                         int _numDir,
                                         int _numContacts) {
  int numOtherConstrs = _A.rows() - _numContacts * (2 + _numDir);
  int n = _numContacts * 3 + numOtherConstrs;
  Eigen::MatrixXd AIntermediate = Eigen::MatrixXd::Zero(n, _A.cols());
  *_AOut = Eigen::MatrixXd::Zero(n, n);
  *_bOut = Eigen::VectorXd::Zero(n);
  int offset = _numDir / 4;
  for (int i = 0; i < _numContacts; ++i) {
    AIntermediate.row(i) = _A.row(i);
    (*_bOut)[i] = _b[i];

    AIntermediate.row(_numContacts + i * 2 + 0) =
        _A.row(_numContacts + i * _numDir + 0);
    AIntermediate.row(_numContacts + i * 2 + 1) =
        _A.row(_numContacts + i * _numDir + offset);
    (*_bOut)[_numContacts + i * 2 + 0] = _b[_numContacts + i * _numDir + 0];
    (*_bOut)[_numContacts + i * 2 + 1] = _b[_numContacts + i * _numDir
                                            + offset];
  }
  for (int i = 0; i < numOtherConstrs; i++) {
    AIntermediate.row(_numContacts * 3 + i) =
        _A.row(_numContacts * (_numDir + 2) + i);
    (*_bOut)[_numContacts * 3 + i] = _b[_numContacts * (_numDir + 2) + i];
  }
  for (int i = 0; i < _numContacts; ++i) {
    _AOut->col(i) = AIntermediate.col(i);
    _AOut->col(_numContacts + i * 2 + 0) =
        AIntermediate.col(_numContacts + i * _numDir + 0);
    _AOut->col(_numContacts + i * 2 + 1) =
        AIntermediate.col(_numContacts + i * _numDir + offset);
  }
  for (int i = 0; i < numOtherConstrs; i++)
    _AOut->col(_numContacts * 3 + i) =
        AIntermediate.col(_numContacts * (_numDir + 2) + i);
}

void LCPSolver::transferSolFromODEFormulation(const Eigen::VectorXd& _x,
                                              Eigen::VectorXd* _xOut,
                                              int _numDir,
                                              int _numContacts) {
  int numOtherConstrs = _x.size() - _numContacts * 3;
  *_xOut = Eigen::VectorXd::Zero(_numContacts * (2 + _numDir)
                                 + numOtherConstrs);

  _xOut->head(_numContacts) = _x.head(_numContacts);

  int offset = _numDir / 4;
  for (int i = 0; i < _numContacts; ++i) {
    (*_xOut)[_numContacts + i * _numDir + 0] = _x[_numContacts + i * 2 + 0];
    (*_xOut)[_numContacts + i * _numDir + offset] = _x[_numContacts + i * 2
                                                       + 1];
  }
  for (int i = 0; i < numOtherConstrs; i++)
    (*_xOut)[_numContacts * (2 + _numDir) + i] = _x[_numContacts * 3 + i];
}

bool LCPSolver::checkIfSolution(const Eigen::MatrixXd& _A,
                                const Eigen::VectorXd& _b,
                                const Eigen::VectorXd& _x) {
  const double threshold = 1e-4;
  int n = _x.size();

  Eigen::VectorXd w = _A * _x + _b;
  for (int i = 0; i < n; ++i) {
    if (w(i) < -threshold || _x(i) < -threshold)
      return false;
    if (abs(w(i) * _x(i)) > threshold)
      return false;
  }
  return true;
}

}  // namespace lcpsolver
}  // namespace dart
