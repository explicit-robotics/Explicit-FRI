/**

The following license terms and conditions apply, unless a redistribution
agreement or other license is obtained by KUKA Deutschland GmbH, Augsburg, Germany.

SCOPE

The software �KUKA Sunrise.Connectivity FRI Client SDK� is targeted to work in
conjunction with the �KUKA Sunrise.Connectivity FastRobotInterface� toolkit.
In the following, the term �software� refers to all material directly
belonging to the provided SDK �Software development kit�, particularly source
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
#ifndef _KUKA_FRI_LBR_JOINT_SINE_OVERLAY_CLIENT_H
#define _KUKA_FRI_LBR_JOINT_SINE_OVERLAY_CLIENT_H

#include "friLBRClient.h"
#include "exp_robots.h"
#include "cart_imp_ctrl.h"
#include "fstream"



/**
 * \brief Test client that can overlay interpolator joint positions with sine waves.
 */
class LBRJointSineOverlayClient : public KUKA::FRI::LBRClient
{
   
public:
      
   /**
    * \brief Constructor.
    * 
    * @param jointMask Bitmask that encodes the joint indices to be overlayed by sine waves
    * @param freqHz Sine frequency in hertz
    * @param amplRad Sine amplitude in radians
    * @param filterCoeff Filter coefficient between 0 (filter off) and 1 (max filter)
    */
   LBRJointSineOverlayClient(unsigned int jointMask, double freqHz, 
         double amplRad, double filterCoeff);
   
   /** 
    * \brief Destructor.
    */
   ~LBRJointSineOverlayClient();
   
   /**
    * \brief Callback for FRI state changes.
    * 
    * @param oldState
    * @param newState
    */
   virtual void onStateChange(KUKA::FRI::ESessionState oldState, KUKA::FRI::ESessionState newState);
   
   /**
    * \brief Callback for the FRI state 'Commanding Active'.
    */
   virtual void command();
      
private:

    // Create iiwa as child of primitive class
    iiwa14 *myLBR;
    
    // Current position and velocity as Eigen vector
    Eigen::VectorXd q;
    Eigen::VectorXd dq;
    
    // External force
    double tau_ext_d[7];
    double qCurr[7];
    Eigen::VectorXd f_ext_0;
    Eigen::VectorXd f_ext_tcp;
    Eigen::VectorXd f_ext_ini;
    Eigen::VectorXd tau_ext;
    
    // DECLARE VARIABLES FOR YOUR CONTROLLER HERE!!!
    Eigen::MatrixXd M;
    Eigen::MatrixXd M_inv;
    Eigen::MatrixXd Lambda_Inv;
    Eigen::MatrixXd Lambda;
    Eigen::MatrixXd H;
    Eigen::MatrixXd J;
    Eigen::MatrixXd J_inv;
    
    // Time parameters for control loop
    double sampleTime;
    double currentTime;
    
    // Files to store data
    std::ofstream File_Fext;
    std::string filenameFile_Fext;
    
   
   int _jointMask;         //!< Bitmask encoding of overlay joints
   double _freqHz;         //!< sine frequency (Hertz)
   double _amplRad;        //!< sine amplitude (radians)
   double _filterCoeff;    //!< filter coefficient
   double _offset;         //!< offset for current interpolation step
   double _phi;            //!< phase of sine wave
   double _stepWidth;      //!< stepwidth for sine 
   
};

#endif // _KUKA_FRI_LBR_JOINT_SINE_OVERLAY_CLIENT_H