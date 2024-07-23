/**

The following license terms and conditions apply, unless a redistribution
agreement or other license is obtained by KUKA Deutschland GmbH, Augsburg, Germany.

SCOPE

The software “KUKA Sunrise.Connectivity FRI Client SDK” is targeted to work in
conjunction with the “KUKA Sunrise.Connectivity FastRobotInterface” toolkit.
In the following, the term “software” refers to all material directly
belonging to the provided SDK “Software development kit”, particularly source
code, libraries, binaries, manuals and technical documentation.

COPYRIGHT

All Rights Reserved
Copyright (C)  2014-2018 
KUKA Deutschland GmbH
Augsburg, Germany

LICENSE 

Redistribution and use of the software in source and binary forms, with or
without modification, are permitted provided that the following conditions are
met:
a) The software is used in conjunction with KUKA products only. 
b) Redistributions of source code must retain the above copyright notice, this
list of conditions and the disclaimer.
c) Redistributions in binary form must reproduce the above copyright notice,
this list of conditions and the disclaimer in the documentation and/or other
materials provided with the distribution. Altered source code of the
redistribution must be made available upon request with the distribution.
d) Modification and contributions to the original software provided by KUKA
must be clearly marked and the authorship must be stated.
e) Neither the name of KUKA nor the trademarks owned by KUKA may be used to
endorse or promote products derived from this software without specific prior
written permission.

DISCLAIMER OF WARRANTY

The Software is provided "AS IS" and "WITH ALL FAULTS," without warranty of
any kind, including without limitation the warranties of merchantability,
fitness for a particular purpose and non-infringement. 
KUKA makes no warranty that the Software is free of defects or is suitable for
any particular purpose. In no event shall KUKA be responsible for loss or
damages arising from the installation or use of the Software, including but
not limited to any indirect, punitive, special, incidental or consequential
damages of any character including, without limitation, damages for loss of
goodwill, work stoppage, computer failure or malfunction, or any and all other
commercial damages or losses. 
The entire risk to the quality and performance of the Software is not borne by
KUKA. Should the Software prove defective, KUKA is not liable for the entire
cost of any service and repair.



\file
\version {1.16}
*/
#include <cstring>
#include <cstdio>

#include <stdio.h>
#include <iostream>
#include <string.h>
#include <math.h>

#include <thread>
#include <vector>
#include <mutex>
#include <unistd.h> // for sleep

#include "LBRJointSineOverlayClient.h"
#include "friLBRState.h"
#include "exp_robots.h"

#include <Python.h>
#include "matplotlibcpp.h"
#include "PlotUtils.h"

std::mutex plotMutex;

using namespace KUKA::FRI;
using namespace std;
namespace plt = matplotlibcpp;

//******************************************************************************
LBRJointSineOverlayClient::LBRJointSineOverlayClient(unsigned int jointMask,
                                                     double freqHz, double amplRad, double filterCoeff)
    : _jointMask(jointMask)
    , _freqHz(freqHz)
    , _amplRad(amplRad)
    , _filterCoeff(filterCoeff)
    , _offset(0.0)
    , _phi(0.0)
    , _stepWidth(0.0)
{
    // Use Explicit-cpp to create your robot
    myLBR = new iiwa14( 1, "Trey");
    myLBR->init( );

    // Current position and velocity
    q  = Eigen::VectorXd::Zero( myLBR->nq );
    dq = Eigen::VectorXd::Zero( myLBR->nq );

    // External force
    for( int i=0; i < myLBR->nq; i++ )
    {
        tau_ext_d[i] = 0.0;
    }

    f_ext_0  = Eigen::VectorXd::Zero( 6 );
    f_ext_tcp = Eigen::VectorXd::Zero( 6 );
    f_ext_ini = Eigen::VectorXd::Zero( 6 );
    tau_ext  = Eigen::VectorXd::Zero( myLBR->nq  );

    // ************************************************************
    // INITIALIZE YOUR VECTORS AND MATRICES HERE
    // ************************************************************
    M = Eigen::MatrixXd::Zero( myLBR->nq, myLBR->nq );
    M_inv = Eigen::MatrixXd::Zero( myLBR->nq, myLBR->nq );
    Lambda_Inv = Eigen::MatrixXd::Zero( 6, 6 );
    Lambda = Eigen::MatrixXd::Zero( 6, 6 );
    H = Eigen::MatrixXd::Zero( 4, 4 );
    J = Eigen::MatrixXd::Zero( 6, myLBR->nq );
    J_inv = Eigen::MatrixXd::Zero( myLBR->nq, 6 );

    // Time variables for control loop
    currentTime = 0;
    sampleTime = 0;

    // File storage
    filenameFile_Fext = "../../prints/File_Fext";
    File_Fext.open(filenameFile_Fext.c_str());

    printf("Ready to rumble!\n");

}


//******************************************************************************
LBRJointSineOverlayClient::~LBRJointSineOverlayClient()
{
    File_Fext.close();
}

//******************************************************************************
void LBRJointSineOverlayClient::onStateChange(ESessionState oldState, ESessionState newState)
{
    LBRClient::onStateChange(oldState, newState);
    // (re)initialize sine parameters when entering Monitoring
    switch (newState)
    {
    case MONITORING_READY:
    {
        _offset = 0.0;
        _phi = 0.0;
        _stepWidth = 2 * M_PI * _freqHz * robotState().getSampleTime();
        sampleTime = robotState().getSampleTime();
        break;
    }
    default:
    {
        break;
    }
    }
}



//******************************************************************************
void LBRJointSineOverlayClient::command()
{

    // ************************************************************
    // get robot measurements
    memcpy( qCurr, robotState().getMeasuredJointPosition(), 7*sizeof(double) );
    memcpy( tau_ext_d, robotState().getExternalTorque(), 7*sizeof(double) );

    for (int i=0; i < myLBR->nq; i++)
    {
        q[i] = qCurr[i];
        tau_ext[i] = tau_ext_d[i];
    }

    // calculate kinematic and dynamic data
    // Homogeneous transformation
    H = myLBR->getForwardKinematics( q );
    Eigen::Vector3d p = H.block< 3, 1 >( 0, 3 );
    Eigen::Matrix3d R = H.block< 3, 3 >( 0, 0 );

    // Jacobian
    J = myLBR->getHybridJacobian( q );

    // mass matrix and its inverse
    M = myLBR->getMassMatrix( q );
    M_inv = this->M.inverse();

    // calculate Lambda least-squares
    double k = 0.01;
//    Lambda_Inv = J * M_inv * J.transpose() + ( k * k ) * Eigen::MatrixXd::Identity( 6, 6 );
    Lambda_Inv = J * M_inv * J.transpose();
    Lambda = Lambda_Inv.inverse();

    // calculate external forces
    J_inv = M_inv * J.transpose() * Lambda;
    f_ext_0 = J_inv.transpose() * tau_ext;

    // reset initial external forces
    if(currentTime < sampleTime)
    {
        f_ext_ini = f_ext_0;
    }
    f_ext_0 = f_ext_0 - f_ext_ini;

    // transform to TCP coordinates
    Eigen::Vector3d f_tcp = Eigen::Vector3d::Zero( 3 );
    Eigen::Vector3d m_tcp = Eigen::Vector3d::Zero( 3 );
    for (int i=0; i < 3; i++)
    {
        f_tcp[i] = f_ext_0[i];
        m_tcp[i] = f_ext_0[i+3];
    }
    f_tcp = R.transpose() * f_tcp;
    m_tcp = R.transpose() * m_tcp;
    for (int i=0; i < 3; i++)
    {
        f_ext_tcp[i] = f_tcp[i];
        f_ext_tcp[i+3] = m_tcp[i];
    }

//    // Print values
//    cout << "************************************************" << endl;
//    cout << "************************************************" << endl;
//    cout << "f_x:";
//    printf("%.2f\n", f_ext_tcp[0]);
//    cout << "f_y:";
//    printf("%.2f\n", f_ext_tcp[1]);
//    cout << "f_z:";
//    printf("%.2f\n", f_ext_tcp[2]);
//    cout << "m_x:";
//    printf("%.2f\n", f_ext_tcp[3]);
//    cout << "m_y:";
//    printf("%.2f\n", f_ext_tcp[4]);
//    cout << "m_z:";
//    printf("%.2f\n", f_ext_tcp[5]);

    // Print values
    cout << "************************************************" << endl;
    cout << "Axial Force:";
    printf("%.2f\n", -f_ext_tcp[2]);

//    // Update shared variables
//    {
//        std::lock_guard<std::mutex> lock(plotMutex);
//        currentTime = this->currentTime;
//        axialForce = f_ext_tcp[2];
//    }

    if( currentTime < sampleTime ){
        File_Fext << "time" << "\t" << "f_x" << "\t" << "f_y" << "\t" << "f_z" << "\t" << "m_x" << "\t" << "m_y" << "\t" << "m_z" << endl;
    }
    File_Fext << currentTime << "\t" << f_ext_tcp[0] << "\t" << f_ext_tcp[1] << "\t" << -f_ext_tcp[2] << "\t" << f_ext_tcp[3] << "\t" << f_ext_tcp[4] << "\t" << f_ext_tcp[5] << endl;

    // no joint offset for now!
    _offset = 0.0;

    // update sample time
    currentTime = currentTime + sampleTime;

    double jointPos[LBRState::NUMBER_OF_JOINTS];
    memcpy(jointPos, robotState().getIpoJointPosition(), LBRState::NUMBER_OF_JOINTS * sizeof(double));
    for (int i=0; i< LBRState::NUMBER_OF_JOINTS; i++)
    {
        if (_jointMask & (1<<i))
        {
            jointPos[i] += _offset;
        }
    }
    robotCommand().setJointPosition(jointPos);

}
//******************************************************************************
// clean up additional defines
#ifdef _USE_MATH_DEFINES
#undef _USE_MATH_DEFINES
#endif
