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
    H = Eigen::Matrix4d::Zero( 4, 4 );
    J = Eigen::MatrixXd::Zero( 6, myLBR->nq );

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

double MyLBRClient::getMaxValue(Eigen::VectorXd myVector)
{
    double auxMax = myVector[0];
    for (int i = 0; i < myVector.size(); i++)
    {
        if (myVector[i] > auxMax)
            auxMax = myVector[i];
    }

    return auxMax;
}


double MyLBRClient::getMinValue(Eigen::VectorXd myVector)
{
    double auxMin = myVector[0];
    for (int i = 0; i < myVector.size(); i++)
    {
        if (myVector[i] < auxMin)
            auxMin = myVector[i];
    }

    return auxMin;
}


Eigen::VectorXd MyLBRClient::addConstraints(Eigen::VectorXd tauStack, double dt)
{
    Eigen::VectorXd dt2 = Eigen::VectorXd::Zero(myLBR->nq,1);
    Eigen::VectorXd dtvar = Eigen::VectorXd::Zero(myLBR->nq,1);
    Eigen::VectorXd qDownBar = Eigen::VectorXd::Zero(myLBR->nq,1);
    Eigen::VectorXd qTopBar = Eigen::VectorXd::Zero(myLBR->nq,1);

    Eigen::VectorXd qDotMaxFromQ = Eigen::VectorXd::Zero(myLBR->nq,1);
    Eigen::VectorXd qDotMinFromQ = Eigen::VectorXd::Zero(myLBR->nq,1);
    Eigen::VectorXd qDotMaxFormQDotDot = Eigen::VectorXd::Zero(myLBR->nq,1);
    Eigen::VectorXd qDotMinFormQDotDot = Eigen::VectorXd::Zero(myLBR->nq,1);
    Eigen::Vector3d vMaxVector = Eigen::Vector3d::Zero(3);
    Eigen::Vector3d vMinVector = Eigen::Vector3d::Zero(3);
    Eigen::VectorXd qDotMaxFinal = Eigen::VectorXd::Zero(myLBR->nq,1);
    Eigen::VectorXd qDotMinFinal = Eigen::VectorXd::Zero(myLBR->nq,1);
    Eigen::VectorXd aMaxqDot = Eigen::VectorXd::Zero(myLBR->nq,1);
    Eigen::VectorXd aMinqDot = Eigen::VectorXd::Zero(myLBR->nq,1);
    Eigen::VectorXd aMaxQ = Eigen::VectorXd::Zero(myLBR->nq,1);
    Eigen::VectorXd aMinQ = Eigen::VectorXd::Zero(myLBR->nq,1);
    Eigen::Vector3d aMaxVector = Eigen::Vector3d::Zero(3);
    Eigen::Vector3d aMinVector = Eigen::Vector3d::Zero(3);
    Eigen::VectorXd qDotDotMaxFinal = Eigen::VectorXd::Zero(myLBR->nq,1);
    Eigen::VectorXd qDotDotMinFinal = Eigen::VectorXd::Zero(myLBR->nq,1);
    Eigen::MatrixXd Iden = Eigen::MatrixXd::Identity(myLBR->nq, myLBR->nq);
    Eigen::VectorXd TauBar = Eigen::VectorXd::Zero(myLBR->nq,1);
    Eigen::VectorXd qDotDotGot = Eigen::VectorXd::Zero(myLBR->nq,1);
    Eigen::MatrixXd Js = Eigen::MatrixXd::Zero(3, myLBR->nq);
    double lowestdtFactor = 10;

    qDownBar = q - myLBR->q_min;
    qTopBar = myLBR->q_max - q;
    dtvar[0] = 3*dt;
    dtvar[1] = 3*dt;
    dtvar[2] = 2*dt;
    dtvar[3] = 3*dt;
    dtvar[4] = dt;
    dtvar[5] = dt;
    dtvar[6] = dt;

    for (int i = 0; i < myLBR->nq; i++)
    {
        //dt2[i] = (lowestdtFactor + (sqrt(lowestdtFactor)*sqrt(10*180/M_PI)))*dt;
        dt2[i] = dtvar[i];
        if (qDownBar[i] < 10*M_PI/180)
        {

            if (qDownBar[i] < 0)
                qDownBar[i] = 0;

            dt2[i] = (lowestdtFactor + (sqrt(lowestdtFactor)*sqrt(qDownBar[i]*180/M_PI)))*dtvar[i];

            if (dt2[i] < lowestdtFactor*dtvar[i])
                dt2[i] = lowestdtFactor*dtvar[i];
        }
        if (qTopBar[i] < 10*M_PI/180)
        {

            if (qTopBar[i] < 0)
                qTopBar[i] = 0;

            dt2[i] = (lowestdtFactor + (sqrt(lowestdtFactor)*sqrt(qTopBar[i]*180/M_PI)))*dtvar[i];
            if (dt2[i] < lowestdtFactor*dtvar[i])
                dt2[i] = lowestdtFactor*dtvar[i];

        }


        qDotMaxFromQ[i] = (myLBR->q_max[i] - q[i]) / dt2[i];
        qDotMinFromQ[i] = (myLBR->q_min[i] - q[i]) / dt2[i];
        qDotMaxFormQDotDot[i] = sqrt(2*myLBR->ddq_max[i]*(myLBR->q_max[i] - q[i]));
        qDotMinFormQDotDot[i] = -sqrt(2*myLBR->ddq_max[i]*(q[i] - myLBR->q_min[i]));

        if(myLBR->q_max[i] - q[i] < 0)
            qDotMaxFormQDotDot[i] = 1000000;

        if(q[i] - myLBR->q_min[i] < 0)
            qDotMinFormQDotDot[i] = -1000000;

        vMaxVector = Eigen::Vector3d(myLBR->dq_max[i], qDotMaxFromQ[i],qDotMaxFormQDotDot[i]);
        qDotMaxFinal[i] = getMinValue(vMaxVector);

        vMinVector = Eigen::Vector3d(myLBR->dq_min[i], qDotMinFromQ[i],qDotMinFormQDotDot[i]);
        qDotMinFinal[i] = getMaxValue(vMinVector);


        aMaxqDot[i] = (qDotMaxFinal[i] - dq[i]) / dtvar[i];
        aMinqDot[i] = (qDotMinFinal[i] - dq[i]) / dtvar[i];
        aMaxQ[i] = 2*(myLBR->q_max[i] - q[i] - dq[i]*dt2[i]) / pow(dt2[i],2);
        aMinQ[i] = 2*(myLBR->q_min[i] - q[i] - dq[i]*dt2[i]) / pow(dt2[i],2);

        aMaxVector = Eigen::Vector3d(aMaxQ[i],aMaxqDot[i],10000000);
        qDotDotMaxFinal[i] = getMinValue(aMaxVector);
        aMinVector = Eigen::Vector3d(aMinQ[i],aMinqDot[i],-10000000);
        qDotDotMinFinal[i] = getMaxValue(aMinVector);

        if  (qDotDotMaxFinal[i] < qDotDotMinFinal[i])
        {
            vMaxVector = Eigen::Vector3d(INFINITY, qDotMaxFromQ[i],qDotMaxFormQDotDot[i]);
            qDotMaxFinal[i] = getMinValue(vMaxVector);

            vMinVector = Eigen::Vector3d(-INFINITY, qDotMinFromQ[i],qDotMinFormQDotDot[i]);
            qDotMinFinal[i] = getMaxValue(vMinVector);

            aMaxqDot[i] = (qDotMaxFinal[i] - dq[i]) / dtvar[i];
            aMinqDot[i] = (qDotMinFinal[i] - dq[i]) / dtvar[i];

            aMaxVector = Eigen::Vector3d(aMaxQ[i],aMaxqDot[i],10000000);
            qDotDotMaxFinal[i] = getMinValue(aMaxVector);
            aMinVector = Eigen::Vector3d(aMinQ[i],aMinqDot[i],-10000000);
            qDotDotMinFinal[i] = getMaxValue(aMinVector);
        }
    }


    Eigen::VectorXd qDotDotS = Eigen::VectorXd::Zero(myLBR->nq);
    Eigen::VectorXd tauS = Eigen::VectorXd::Zero(myLBR->nq);
    Eigen::MatrixXd Psat = Iden;
    bool LimitedExceeded = true;
    bool CreateTaskSat = false;
    int NumSatJoints = 0;
    Eigen::VectorXd theMostCriticalOld = Eigen::VectorXd::Zero(myLBR->nq);
    theMostCriticalOld.conservativeResize(1);
    theMostCriticalOld[0] = 100;
    bool isThere = false;
    int iO = 0;
    int cycle = 0;
    while (LimitedExceeded == true)
    {

        LimitedExceeded = false;

        if (CreateTaskSat == true)
        {
            Js.conservativeResize(NumSatJoints,myLBR->nq);
            for (int i = 0; i < NumSatJoints; i++)
            {
                for(int k = 0;k < myLBR->nq; k++)
                {
                    Js(i,k) = 0;
                }
                int m = theMostCriticalOld[i];
                Js(i,m) = 1;
            }

            Eigen::MatrixXd Minv = M.inverse();
            Eigen::MatrixXd LambdaSatInv = Js * Minv * Js.transpose();
            Eigen::MatrixXd LambdaSatInv_aux = LambdaSatInv*LambdaSatInv.transpose();
            Eigen::MatrixXd LambdaSat_aux = LambdaSatInv_aux.inverse();
            Eigen::MatrixXd LambdaSat = LambdaSatInv.transpose()*LambdaSat_aux;

            Eigen::MatrixXd JsatBar = Minv * Js.transpose()*LambdaSat;
            Psat = Iden - Js.transpose()*JsatBar.transpose();
            Eigen::VectorXd xDotDot_s = Js*qDotDotS; //+ JsDot*myRobot->qDot;
            tauS = Js.transpose()*(LambdaSat*xDotDot_s);


        }

        TauBar = tauS + Psat*tauStack;
        qDotDotGot = M.inverse() * TauBar; // it should -g -c

        isThere = false;
        for (int i = 0; i < myLBR->nq; i++)
        {
            if ((qDotDotMaxFinal[i] + 0.001 < qDotDotGot[i])  || (qDotDotGot[i] < qDotDotMinFinal[i] - 0.001))
            {
                LimitedExceeded = true;
                CreateTaskSat = true;

                for (int k = 0; k < theMostCriticalOld.size(); k++)
                {
                    if (i == theMostCriticalOld[k])
                    {
                        isThere = true;
                    }
                }
                if (isThere == false)
                {

                    theMostCriticalOld.conservativeResize(iO + 1);
                    theMostCriticalOld[iO] = i;
                    iO += 1;
                }
            }
        }

        if (LimitedExceeded == true)
        {
            NumSatJoints = iO;
            theMostCriticalOld.conservativeResize(iO);
            cycle +=1;
            if (cycle > 8)
                LimitedExceeded = false;

            for (int i = 0; i < theMostCriticalOld.size(); i++)
            {
                int jM = theMostCriticalOld[i];

                if (qDotDotGot[jM] > qDotDotMaxFinal[jM])
                    qDotDotS[jM] = qDotDotMaxFinal[jM];

                if (qDotDotGot[jM] < qDotDotMinFinal[jM])
                    qDotDotS[jM] = qDotDotMinFinal[jM];
            }
        }
    }

    Eigen::VectorXd SJSTorque = TauBar;
    return SJSTorque;
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

    for (int i=0; i < myLBR->nq; i++)
    {
        q[i] = qCurr[i];
    }

    for (int i=0; i < 7; i++)
    {
        dq[i] = (qCurr[i] - qOld[i]) / sampleTime;
    }

    // ************************************************************
    // Calculate kinematics and dynamics

    //if(currentTime < sampleTime)
    //{
    // To get initial values for first time you are in the control loop
    //}

    M = myLBR->getMassMatrix( q );
    M( 6, 6 ) = 45 * M( 6, 6 );

    // ************************************************************
    // YOUR CONTROLLER COMES HERE!
    // ************************************************************


    // ************************************************************
    // Control torque

    // tau_motion = ...
    // Include joint limits
    tau_motion = Eigen::VectorXd::Zero( myLBR->nq );
    tau_motion = addConstraints(tau_motion,0.004);

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
