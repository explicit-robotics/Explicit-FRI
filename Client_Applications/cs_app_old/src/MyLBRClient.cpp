/**

DISCLAIMER OF WARRANTY

The Software is provided "AS IS" and "WITH ALL FAULTS,"
without warranty of any kind, including without limitation the warranties
of merchantability, fitness for a particular purpose and non-infringement.
KUKA makes no warranty that the Software is free of defects or is suitable
for any particular purpose. In no event shall KUKA be responsible for loss
or damages arising from the installation or use of the Software,
including but not limited to any indirect, punitive, special, incidental
or consequential damages of any character including, without limitation,
damages for loss of goodwill, work stoppage, computer failure or malfunction,
or any and all other commercial damages or losses.
The entire risk to the quality and performance of the Software is not borne by KUKA.
Should the Software prove defective, KUKA is not liable for the entire cost
of any service and repair.


COPYRIGHT

All Rights Reserved
Copyright (C)  2014-2015
KUKA Roboter GmbH
Augsburg, Germany

This material is the exclusive property of KUKA Roboter GmbH and must be returned
to KUKA Roboter GmbH immediately upon request.
This material and the information illustrated or contained herein may not be used,
reproduced, stored in a retrieval system, or transmitted in whole
or in part in any way - electronic, mechanical, photocopying, recording,
or otherwise, without the prior written consent of KUKA Roboter GmbH.






\file
\version {1.9}
*/
#include <stdio.h>
#include <iostream>
#include <string.h>
#include <math.h>
#include <chrono>

#include "MyLBRClient.h"
#include "exp_robots.h"

#include <boost/thread.hpp>
#include <boost/chrono.hpp>

#include "my_diff_jacobians.h"

using namespace std;

#ifndef M_PI
#define M_PI 3.14159265358979
#endif

#ifndef NCoef
#define NCoef 1
#endif

static double filterOutput[7][NCoef+1]; //output samples. Static variables are initialised to 0 by default.
static double filterInput[7][NCoef+1]; //input samples. Static variables are initialised to 0 by default.


//******************************************************************************
MyLBRClient::MyLBRClient(double freqHz, double amplitude){

    /** Initialization */

    // THIS CONFIGURATION MUST BE THE SAME AS FOR THE JAVA APPLICATION!!
    qInitial[0] = 32.54 *M_PI/180;
    qInitial[1] = 56.74 *M_PI/180;
    qInitial[2] = 0.0   *M_PI/180;
    qInitial[3] = -65.98 *M_PI/180;
    qInitial[4] = 0.0   *M_PI/180;
    qInitial[5] = 57.29 *M_PI/180;
    qInitial[6] = 32.55  *M_PI/180;

    // Use Explicit-cpp to create your robot
    myLBR = new iiwa14( 1, "Trey");
    myLBR->init( );

    // Current position and velocitz
    q  = Eigen::VectorXd::Zero( myLBR->nq );
    dq = Eigen::VectorXd::Zero( myLBR->nq );

    // Time variables for control loop
    currentTime = 0;
    sampleTime = 0;

    // Initialize joint torques and joint positions (also needed for waitForCommand()!)
    for( int i=0; i < myLBR->nq; i++ )
    {
        qCurr[i] = qInitial[i];
        qOld[i] = qInitial[i];
        qApplied[i] = 0.0;
        torques[i] = 0.0;
    }

    tau_motion    = Eigen::VectorXd::Zero( myLBR->nq );
    tau_previous  = Eigen::VectorXd::Zero( myLBR->nq );
    tau_prev_prev = Eigen::VectorXd::Zero( myLBR->nq );
    tau_total     = Eigen::VectorXd::Zero( myLBR->nq );


    // ************************************************************
    // INITIALIZE YOUR VECTORS AND MATRICES HERE
    // ************************************************************
    M = Eigen::MatrixXd::Zero( myLBR->nq, myLBR->nq );
    M_inv = Eigen::MatrixXd::Zero( myLBR->nq, myLBR->nq );
    H = Eigen::Matrix4d::Zero( 4, 4 );
    J = Eigen::MatrixXd::Zero( 6, myLBR->nq );
    J_inv = Eigen::MatrixXd::Zero( myLBR->nq, 6 );

    tauExt  = Eigen::VectorXd::Zero( myLBR->nq  );
    f_ext_0  = Eigen::VectorXd::Zero( 6 );

    Kp = Eigen::MatrixXd::Identity( 3, 3 );
    Kp = 400 * Kp;
    K = Eigen::MatrixXd::Zero( 6, 6 );
    for (int i=0; i<3; i++)
    {
        K(i,i) = Kp(i,i);
    }

    Kr = Eigen::MatrixXd::Identity( 3, 3 );
    Kr = 20 * Kr;
    for (int i=0; i<3; i++)
    {
        K(i+3,i+3) = Kr(i,i);
    }

    Kq = Eigen::MatrixXd::Identity( 7, 7 );
    Kq = 10 * Kq;
    Bq = Eigen::MatrixXd::Identity( 7, 7 );
    Bq = 5 * Bq;
    Kq_thread = Eigen::MatrixXd::Zero( 7, 7 );
    Bq_thread = Eigen::MatrixXd::Zero( 7, 7 );
    K_fin_thread = Eigen::MatrixXd::Zero( 7, 7 );
    Kq_thread = Eigen::MatrixXd::Zero( 7, 7 );
    q_0 = Eigen::VectorXd::Zero( myLBR->nq );

    mutex.unlock();

    // Initialize Impedance Parameters
    myImpedanceParams = new ImpedanceParameters();
    myImpedanceParams->K_q = Kq;
    myImpedanceParams->B_q = Bq;

    // Initial print
    printf( "Exp[licit](c)-cpp-FRI, https://explicit-robotics.github.io \n\n" );
    printf( "Robot '" );
    printf( "%s", myLBR->Name );
    printf( "' initialised. Ready to rumble! \n\n" );

}


/**
* \brief Destructor
*
*/
MyLBRClient::~MyLBRClient()
{


}


/**
* \brief Implements an IIR Filter which is used to send the previous joint position to the command function, so that KUKA's internal friction compensation can be activated. The filter function was generated by the application WinFilter (http://www.winfilter.20m.com/).
*
* @param NewSample The current joint position to be provided as input to the filter.
*/
void iir(double NewSample[7])
{
    double ACoef[ NCoef+1 ] = {
        0.05921059165970496400,
        0.05921059165970496400
    };

    double BCoef[ NCoef+1 ] = {
        1.00000000000000000000,
        -0.88161859236318907000
    };

    int n;

    // Shift the old samples
    for ( int i=0; i<7; i++ )
    {
        for( n=NCoef; n>0; n-- )
        {
            filterInput[i][n] = filterInput[i][n-1];
            filterOutput[i][n] = filterOutput[i][n-1];
        }
    }

    // Calculate the new output
    for (int i=0; i<7; i++)
    {
        filterInput[i][0] = NewSample[i];
        filterOutput[i][0] = ACoef[0] * filterInput[i][0];
    }

    for (int i=0; i<7; i++)
    {
        for(n=1; n<=NCoef; n++)
            filterOutput[i][0] += ACoef[n] * filterInput[i][n] - BCoef[n] * filterOutput[i][n];
    }
}

//******************************************************************************
void MyLBRClient::onStateChange(ESessionState oldState, ESessionState newState)
{
    LBRClient::onStateChange(oldState, newState);
    // react on state change events
    switch (newState)
    {
    case MONITORING_WAIT:
    {
        break;
    }
    case MONITORING_READY:
    {
        sampleTime = robotState().getSampleTime();
        break;
    }
    case COMMANDING_WAIT:
    {
        break;
    }
    case COMMANDING_ACTIVE:
    {
        break;
    }
    default:
    {
        break;
    }
    }
}

//******************************************************************************
void MyLBRClient::monitor()
{

    // Copied from FRIClient.cpp
    robotCommand().setJointPosition(robotState().getCommandedJointPosition());

    // Copy measured joint positions (radians) to _qcurr, which is a double

    memcpy( qCurr, robotState().getMeasuredJointPosition(), 7*sizeof(double) );

    // Initialise the q for the previous NCoef timesteps

    for( int i=0; i<NCoef+1; i++ )
    {
        iir(qCurr);
    }
}

//******************************************************************************
void MyLBRClient::waitForCommand()
{
    // If we want to command torques, we have to command them all the time; even in
    // waitForCommand(). This has to be done due to consistency checks. In this state it is
    // only necessary, that some torque vlaues are sent. The LBR does not take the
    // specific value into account.

    if(robotState().getClientCommandMode() == TORQUE){

        robotCommand().setTorque(torques);
        robotCommand().setJointPosition(robotState().getIpoJointPosition());            // Just overlaying same position
    }

}

//******************************************************************************
void MyLBRClient::command()
{

    // ************************************************************
    // Get robot measurements

    memcpy( qOld, qCurr, 7*sizeof(double) );
    memcpy( qCurr, robotState().getMeasuredJointPosition(), 7*sizeof(double) );
    memcpy( tauExternal, robotState().getExternalTorque(), 7*sizeof(double) );

    for (int i=0; i < myLBR->nq; i++)
    {
        q[i] = qCurr[i];
    }

    for (int i=0; i < 7; i++)
    {
        dq[i] = (qCurr[i] - qOld[i]) / sampleTime;
    }

    for (int i=0; i < 7; i++)
    {
        tauExt[i] = tauExternal[i];
    }


    // ************************************************************
    // Calculate kinematics and dynamics

    if(currentTime < sampleTime)
    {
        // To get initial values for first time you are in the control loop
        q_0 = q;
    }

    // Jacobian
    J = myLBR->getHybridJacobian( q );

    // Mass matrix
    M = myLBR->getMassMatrix( q );
    M( 6, 6 ) = 40 * M( 6, 6 );
    M_inv = this->M.inverse();

    // Lambda least-squares
    double k = 0.01;
    Eigen::MatrixXd Lambda_Inv = J * M_inv * J.transpose() + ( k * k ) * Eigen::MatrixXd::Identity( 6, 6 );
    //Eigen::MatrixXd Lambda_Inv = J * M_inv * J.transpose();
    Eigen::MatrixXd Lambda = Lambda_Inv.inverse();

    // External forces
    J_inv = M_inv * J.transpose() * Lambda;
    f_ext_0 = J_inv.transpose() * tauExt;


    // ************************************************************
    // Impedance optimization

    // INITIALIZE K AND B BEFORE FIRST THREAD IS FINISHED!

    if(currentTime < sampleTime)
    {
        impedanceOptimizerThread = boost::thread(&MyLBRClient::impedanceOptimizer, this);
        impedanceOptimizerThread.detach();
    }

    mutex.lock();
    // CHECK THIS!!!!!!!!!!!!!!!!!!!!
    Kq_thread = Kq_thread;
    Bq_thread = Bq_thread;
    K_kin_thread = K_kin_thread;

    mutex.unlock();

    Kq = Kq_thread;
    Bq = Bq_thread;
    Eigen::MatrixXd K_kin = K_kin_thread;

    Eigen::MatrixXd Kasym = Kq - Kq.transpose( );
    Eigen::MatrixXd Basym = Bq - Bq.transpose( );

    //    cout << "K_kin:" << endl;
    //    cout << K_kin << endl;
    //    cout << "Kq:" << endl;
    //    cout << Kq << endl;
    //    cout << "Bq:" << endl;
    //    cout << Bq << endl;

    //    cout << "J: " << endl;
    //    cout << J << endl;

    //    cout << "J^TKJ: " << endl;
    //    cout << J.transpose() * K * J << endl;

    //    cout << "K asymmetric part:" << endl;
    //    cout << Kasym.maxCoeff( ) << endl;
    //    cout << "B asymmetric part:" << endl;
    //    cout << Basym.maxCoeff( ) << endl;

    myImpedanceParams->K_q = Kq;
    myImpedanceParams->B_q = Bq;


    // ************************************************************
    // Control torque
    tau_motion = myImpedanceParams->K_q * (q_0 - q) - myImpedanceParams->B_q * dq;

    //    cout << "tau_motion:" << endl;
    //    cout << tau_motion << endl;

    // Just gravity compensation
    //tau_motion = Eigen::VectorXd::Zero( myLBR->nq );

    // Include joint limits
    //tau_motion = myLBR->addIIWALimits( q, dq, M, tau_motion, 0.004 );


    // ************************************************************
    // YOUR CODE ENDS HERE!
    // ************************************************************

    // A simple filter for the torque command
    tau_total = ( tau_motion + tau_previous + tau_prev_prev ) / 3;

    for ( int i=0; i<7; i++ )
    {
        qApplied[i] = filterOutput[i][0];
        torques[i] = tau_total[i];
    }

    // Command values (must be double arrays!)
    if (robotState().getClientCommandMode() == TORQUE)
    {
        robotCommand().setJointPosition(qApplied);
        robotCommand().setTorque(torques);
    }

    // IIR filter input
    iir(qCurr);

    // Update
    if (currentTime == 0.0)
    {
        tau_previous = tau_motion;
        tau_prev_prev = tau_motion;
    }
    tau_previous = tau_motion;
    tau_prev_prev = tau_previous;

    currentTime = currentTime + sampleTime;

}

//******************************************************************************
void MyLBRClient::impedanceOptimizer()
{

    while(true){

        // Calculate standard Hessian
        Eigen::MatrixXd kq = J.transpose() * K * J;

        cout << "kq:" << endl;
        cout << kq << endl;

        // Christoffel symbols based on Kinematic Connection for Hybrid Jacobian
        Eigen::MatrixXd CS = Eigen::MatrixXd::Zero( 6, 6 );
        CS( 3, 4 ) = f_ext_0( 5 )/2;
        CS( 3, 5 ) = -f_ext_0( 4 )/2;
        CS( 4, 3 ) = -f_ext_0( 5 )/2;
        CS( 4, 5 ) = f_ext_0( 3 )/2;
        CS( 5, 3 ) = f_ext_0( 4 )/2;
        CS( 5, 4 ) = -f_ext_0( 3 )/2;

        // Calculate Jacobian Partials
        Eigen::MatrixXd dJH1 = dJH_T_dq1( q[ 0 ], q[ 1 ], q[ 2 ], q[ 3 ], q[ 4 ], q[ 5 ], q[ 6 ] );
        Eigen::MatrixXd dJH2 = dJH_T_dq2( q[ 0 ], q[ 1 ], q[ 2 ], q[ 3 ], q[ 4 ], q[ 5 ], q[ 6 ] );
        Eigen::MatrixXd dJH3 = dJH_T_dq3( q[ 0 ], q[ 1 ], q[ 2 ], q[ 3 ], q[ 4 ], q[ 5 ], q[ 6 ] );
        Eigen::MatrixXd dJH4 = dJH_T_dq4( q[ 0 ], q[ 1 ], q[ 2 ], q[ 3 ], q[ 4 ], q[ 5 ], q[ 6 ] );
        Eigen::MatrixXd dJH5 = dJH_T_dq5( q[ 0 ], q[ 1 ], q[ 2 ], q[ 3 ], q[ 4 ], q[ 5 ], q[ 6 ] );
        Eigen::MatrixXd dJH6 = dJH_T_dq6( q[ 0 ], q[ 1 ], q[ 2 ], q[ 3 ], q[ 4 ], q[ 5 ], q[ 6 ] );
        Eigen::MatrixXd dJH7 = dJH_T_dq7( q[ 0 ], q[ 1 ], q[ 2 ], q[ 3 ], q[ 4 ], q[ 5 ], q[ 6 ] );

        // Calculate Kinematic Stiffness
        Eigen::MatrixXd K_kin = Eigen::MatrixXd::Zero( 7, 7 );

        K_kin.col( 0 ) = dJH1 * f_ext_0;
        K_kin.col( 1 ) = dJH2 * f_ext_0;
        K_kin.col( 2 ) = dJH3 * f_ext_0;
        K_kin.col( 3 ) = dJH4 * f_ext_0;
        K_kin.col( 4 ) = dJH5 * f_ext_0;
        K_kin.col( 5 ) = dJH6 * f_ext_0;
        K_kin.col( 6 ) = dJH7 * f_ext_0;

        K_kin = K_kin + J.transpose() * CS * J;
        Eigen::MatrixXd K_final = K_kin + kq;

        // Calculate Damping
        Eigen::MatrixXd sqrt_M = Eigen::MatrixXd::Zero( 7, 7 );
        for(int i=0; i<7; i++)
        {
            for(int j=0; j<7; j++)
            {
                if( M( i, j ) < 0 )
                {
                    M( i, j ) = 0.0;
                }
                sqrt_M( i, j ) = std::sqrt( M( i, j ) );

                if( std::isnan( sqrt_M( i, j ) ) )
                {
                    sqrt_M( i, j ) = 0.0;
                }

            }
        }

        Eigen::MatrixXd sqrt_K = Eigen::MatrixXd::Zero( 7, 7 );
        for(int i=0; i<7; i++)
        {
            for(int j=0; j<7; j++)
            {
                if( kq( i, j ) < 0 )
                {
                    kq( i, j ) = 0.0;
                }

                sqrt_K( i, j ) = std::sqrt( kq( i, j ) );

                if( std::isnan( sqrt_K( i, j ) ) )
                {
                    sqrt_K( i, j ) = 0.0;
                }
            }
        }

        Eigen::MatrixXd zeta = Eigen::MatrixXd::Identity( 7, 7 );
        zeta = 0.7 * zeta;
        Eigen::MatrixXd bq = sqrt_M * zeta * sqrt_K + sqrt_K * zeta * sqrt_M;


        //****************** Update everyting at the end with one Mutex ******************//
        mutex.lock();

        Kq_thread = K_final;
        //        Kq_thread = kq;
        Bq_thread = bq;
        K_kin_thread = K_kin;

        mutex.unlock();

    }

}
