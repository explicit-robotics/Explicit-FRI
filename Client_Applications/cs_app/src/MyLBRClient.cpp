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
// A Function to read a CSV file and save it as a 2D matrix
// Dimensions for row x col are 7 x N, where 7 (row) is the degrees of freedom of the KUKA, N (col) is the number of data points
Eigen::MatrixXd readCSV( const string& filename )
{
    ifstream file( filename );
    if( !file.is_open( ) )
    {
        cerr << "Error: Couldn't open the file: " << filename << endl;
        exit( 1 );
    }

    // Read the values as 2D vector array
    vector<vector<double>> values;

    string line;
    int lineNum = 0;
    int numCols = 0;

    while ( getline( file, line ) )
    {
        ++lineNum;
        stringstream ss( line );
        string cell;
        vector<double> row;
        while ( getline( ss, cell, ',' ) )
        {
            try
            {
                row.push_back( stod( cell ) );
            } catch ( const std::invalid_argument& e )
            {
                cerr << "Error: Invalid argument at line " << lineNum << ", column: " << row.size() + 1 << endl;
                exit( 1 );
            }
        }

        values.push_back( row );
        if ( numCols == 0 )
        {
            numCols = row.size( );
        }
        else if ( row.size() != numCols )
        {
            cerr << "Error: Inconsistent number of columns in the CSV file." << endl;
            exit( 1 );
        }
    }

    // Error if CSV File is empty
    if ( values.empty( ) )
    {
        cerr << "Error: CSV file is empty." << endl;
        exit(1);
    }

    // Create Eigen Matrix
    Eigen::MatrixXd mat( values.size(), numCols );
    for (int i = 0; i < values.size(); i++)
    {
        for (int j = 0; j < numCols; j++)
        {
            mat(i, j) = values[i][j];
        }
    }
    return mat;
}


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

    // Current joint configuration and velocity
    q  = Eigen::VectorXd::Zero( myLBR->nq );
    dq = Eigen::VectorXd::Zero( myLBR->nq );

    // Desired joint configuration
    q_0 = Eigen::VectorXd::Zero( myLBR->nq );

    // Time variables for control loop
    currentTime = 0;
    sampleTime = 0;
    n_step = 0;         // Number of current time steps
    t_pressed = 0;      // Time after the button push

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

    // Read the joint value data
    // FILL-IN YOUR CSV FILE ADDRESS!!
    q_data = readCSV( "/home/newman_lab/Desktop/Explicit-FRI/Client_Applications/cs_app/data/kin_stiff_test1.csv" );

    // Number of joint-trajectory data points, and its current iteration
    N_data = q_data.cols( );
    N_curr = 0;

    // Get the initial joint position from the data
    q_init1 = q_data.col( 0 );

    // Also get the initial joint position of current joint position
    q_init0 = Eigen::VectorXd::Zero( myLBR->nq );
    for( int i = 0; i < myLBR->nq; i++ )
    {
        q_init0( i ) = qInitial[ i ];
    }

    // Generate a minimum-jerk trajectory in joint-space to move to the starting location
    mjt_q = new MinimumJerkTrajectory( 7, q_init0, q_init1, 5.0, 1.0 );

    // Boolean variable of button pressed
    is_pressed = false;


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

    // Cartesian stiffness
    Kp = Eigen::MatrixXd::Identity( 3, 3 );
    Kp = 850 * Kp;
    K = Eigen::MatrixXd::Zero( 6, 6 );
    for ( int i=0; i<3; i++ )
    {
        K( i, i ) = Kp( i, i );
    }

    Kr = Eigen::MatrixXd::Identity( 3, 3 );
    Kr = 50 * Kr;
    for (int i=0; i<3; i++)
    {
        K( i+3, i+3 ) = Kr( i, i );
    }

    // Joint impedances
    Kq = Eigen::MatrixXd::Identity( 7, 7 );
    Kq = 100 * Kq;
    Bq = Eigen::MatrixXd::Identity( 7, 7 );
    Bq = 10 * Bq;
    Kq_last = Kq;
    Bq_last = Bq;

    Kq_thread = Kq;
    Bq_thread = Bq;

    K_kin = Eigen::MatrixXd::Zero(7,7);
    K_fin_thread = Eigen::MatrixXd::Zero( 7, 7 );
    Kq_thread = Eigen::MatrixXd::Zero( 7, 7 );

    // Lock mutex initially
    mutex.unlock();

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

    // Calculate the virtual trajectory
    q_0 = mjt_q->getPosition( currentTime );

    // Add the delta joint q0 position via data, to replay the variable
    // If button pressed, and after 2-seconds
    if( is_pressed && t_pressed >= 2 )
    {
        q_0 += q_data.col( N_curr ) - q_init1;
    }


    // ************************************************************
    // Calculate kinematics and dynamics

    // Jacobian
    J = myLBR->getHybridJacobian( q );

    // Mass matrix
    M = myLBR->getMassMatrix( q );
    M( 6, 6 ) = 40 * M( 6, 6 );
    M_inv = this->M.inverse();

    if(currentTime < sampleTime)
    {
        // Joint stiffness from initial configuration to start of recorded joint path
        Kq = J.transpose() * K * J;

        // Handle nullspace through simple least squares method
        Kq = Kq + 20 * Eigen::MatrixXd::Identity(7,7);
    }

    // ************************************************************
    // Calculate external forces

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

    if(currentTime < sampleTime)
    {
        impedanceOptimizerThread = boost::thread(&MyLBRClient::impedanceOptimizer, this);
        impedanceOptimizerThread.detach();
    }

    // Store last impedance values
    Kq_last = Kq;
    Bq_last = Bq;

    // Only update impedances for recorded trajectory
    if( is_pressed && t_pressed >= 2 )
    {
        mutex.lock();

        Kq = Kq_thread;
        Bq = Bq_thread;
        K_kin = K_kin_thread;

        mutex.unlock();
    }

    // Safety check (DELETE LATER!)
    for (int i=0; i < 7; i++)
    {
        for (int j=0; j < 7; j++)
        {
            if( Kq(i,j) >= 1e+5 )
            {
                Kq(i,j) = Kq_last(i,j);
                cout << "************************************************HERE*************************************" << endl;
                cout << "************************************************HERE*************************************" << endl;
                cout << "************************************************HERE*************************************" << endl;
                cout << "************************************************HERE*************************************" << endl;
                cout << "************************************************HERE*************************************" << endl;
            }

            if( Bq(i,j) >= 1e+5 )
            {
                Bq(i,j) = Bq_last(i,j);
            }
        }
    }

    // Just some plots for now
    Eigen::MatrixXd Kasym = Kq - Kq.transpose( );
    Eigen::MatrixXd Basym = Bq - Bq.transpose( );

    //    cout << "K_kin:" << endl;
    //    cout << K_kin << endl;
    cout << "Kq:" << endl;
    cout << Kq << endl;
    cout << "Bq:" << endl;
    cout << Bq << endl;

    //    cout << "J: " << endl;
    //    cout << J << endl;

    //    cout << "J^TKJ: " << endl;
    //    cout << J.transpose() * K * J << endl;

    //    cout << "K asymmetric part:" << endl;
    //    cout << Kasym.maxCoeff( ) << endl;
    //    cout << "B asymmetric part:" << endl;
    //    cout << Basym.maxCoeff( ) << endl;


    // ************************************************************
    // Control torque
    tau_motion = Kq * (q_0 - q) - Bq * dq;

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

    // Check if Button Pressed
    if( robotState().getBooleanIOValue( "MediaFlange.UserButton" ) && !is_pressed )
    {
        is_pressed = true;

        // Turn on imitation learning
        cout << "Button Pressed!" << std::endl;
    }

    // If button pressed
    if( is_pressed )
    {
        // Update time for t_pressed
        t_pressed += sampleTime;

        // Change 1 to 2 to make it slower
        if ( t_pressed >= 2 )
        {
            // Update rate of n_step.
            if (  ( n_step % 1 ) == 0 )
            {
                N_curr++;

                // The N_curr must not exceed the column size.
                if( N_curr >= N_data )
                {
                    N_curr = N_data-1;
                }
            }
        }

    }

    n_step++;
    currentTime = currentTime + sampleTime;

}

//******************************************************************************
void MyLBRClient::impedanceOptimizer()
{

    while(true){

        // Calculate standard Hessian
        Eigen::MatrixXd kq = J.transpose() * K * J;

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

        // Damping design
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
        Bq_thread = bq;
        K_kin_thread = K_kin;

        mutex.unlock();

    }

}
