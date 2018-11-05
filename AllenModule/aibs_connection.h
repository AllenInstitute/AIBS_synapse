/*
 *  drop_odd_spike_connection.h
 *
 *  This file is part of NEST.
 *
 *  Copyright (C) 2004 The NEST Initiative
 *
 *  NEST is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  NEST is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with NEST.  If not, see <http://www.gnu.org/licenses/>.
 *
 */

#ifndef AIBS_CONNECTION_H
#define AIBS_CONNECTION_H

// Includes from nestkernel:
#include "connection.h"
#include <iostream>

/* BeginDocumentation
  Name: drop_odd_spike - Synapse dropping spikes with odd time stamps.

  Description:
  This synapse will not deliver any spikes with odd time stamps, while spikes
  with even time stamps go through unchanged.

  Transmits: SpikeEvent

  Remarks:
  This synapse type is provided only for illustration purposes in MyModule.

  SeeAlso: synapsedict
*/

namespace allennest
{
inline double dndt(double t, double n, double Tau_r)
{
    if (Tau_r>1e-10)
    {
      return (1-n)/Tau_r;
    }
    else
    {
      return 0;
    }
}

inline double dpdt(double t, double p, double p0, double Tau_f)
{
    if (Tau_f>1e-10)
    {
      return (p0-p)/Tau_f;
    }
    else
    {
      return 0;
    }
}

inline double dtrdt(double t, double Tau_r, double Tau_FDR, double Tau_r0)
{
     if (Tau_FDR>1e-10)
     {
       return (Tau_r0-Tau_r)/Tau_FDR;
     }
     else
     {
       return 0;
     }

}
inline double dSdt(double t, double S, double Tau_D)
{
    if (Tau_D>1e-10)
    {    
      return (1-S)/Tau_D;
    }
    else
    {
      return 0;
    }
}

inline double dp0dt(double t, double p0, double p0bar, double Tau_i)
{
    if (Tau_i>1e-10)
    {
      return (p0bar-p0)/Tau_i;
    }
    else
    {
      return 0;
    }
}

/**
 * Connection class for illustration purposes.
 *
 * For a discussion of how synapses are created and represented in NEST 2.6,
 * please see Kunkel et al, Front Neuroinform 8:78 (2014), Sec 3.3.
 */

// This file relies on inline functions (such as dn/dt) in ai_connection.h. If this file is installed alone, please copy and paste them here.
template < typename targetidentifierT >
class AibsConnection : public nest::Connection< targetidentifierT >
{
private:
  double weight_;
  double tau_f_; //!< [ms] time constant for fascilitation
  double tau_r_; //!< [ms] time constant for recovery
  double tau_fdr_; //!< [ms] time constant for baseline of tau_rec
  double tau_r0_; //!< [ms] asymptotic value of baseline of tau_rec
  double p0_;       //!< asymptotic value of probability of release
  double n_;       //!< amount of resources in recovered state
  double p_;       //!< actual probability of release
  double a_fdr_;       //!< fraction of changes in tau_rec
  double S_;       //the fraction of non-desensitized receptors
  double tau_d_;   //the 
  double a_D_;
  double a_i_;
  double p0bar_;
  double tau_i_;
  double flag_f;
  double flag_Tau_r;
  double flag_S;
  double flag_p0;
  double t_lastspike_; //!< time point of last spike emitted
public:
  //! Type to use for representing common synapse properties
  typedef nest::CommonSynapseProperties CommonPropertiesType;

  //! Shortcut for base class
  typedef nest::Connection< targetidentifierT > ConnectionBase;

  /**
   * Default Constructor.
   * Sets default values for all parameters. Needed by GenericConnectorModel.
   */
  AibsConnection();  
  
  AibsConnection(const AibsConnection&);  
  
  //! Default Destructor.
  ~AibsConnection()
  {
  }
  // Explicitly declare all methods inherited from the dependent base
  // ConnectionBase. This avoids explicit name prefixes in all places these
  // functions are used. Since ConnectionBase depends on the template parameter,
  // they are not automatically found in the base class.
  using ConnectionBase::get_delay_steps;
  using ConnectionBase::get_delay;
  using ConnectionBase::get_rport;
  using ConnectionBase::get_target;


  /**
   * Helper class defining which types of events can be transmitted.
   *
   * These methods are only used to test whether a certain type of connection
   * can be created.
   *
   * `handles_test_event()` should be added for all event types that the
   * synapse can transmit. The methods shall return `invalid_port_`; the
   * return value will be ignored.
   *
   * Since this is a synapse model dropping spikes, it is only for spikes,
   * therefore we only implement `handles_test_event()` only for spike
   * events.
   *
   * See Kunkel et al (2014), Sec 3.3.1, for background information.
   */
  class ConnTestDummyNode : public nest::ConnTestDummyNodeBase
  {
  public:
    using nest::ConnTestDummyNodeBase::handles_test_event;
    nest::port
    handles_test_event( nest::SpikeEvent&, nest::rport )
    {
      return nest::invalid_port_;
    }

    nest::port
    handles_test_event( nest::DSSpikeEvent&, nest::rport )
    {
      return nest::invalid_port_;
    }
  };

  /**
   * Check that requested connection can be created.
   *
   * This function is a boilerplate function that should be included unchanged
   * in all synapse models. It is called before a connection is added to check
   * that the connection is legal. It is a wrapper that allows us to call
   * the "real" `check_connection_()` method with the `ConnTestDummyNode
   * dummy_target;` class for this connection type. This avoids a virtual
   * function call for better performance.
   *
   * @param s  Source node for connection
   * @param t  Target node for connection
   * @param receptor_type  Receptor type for connection
   * @param lastspike Time of most recent spike of presynaptic (sender) neuron,
   *                  not used here
   */
  void
  check_connection( nest::Node& s,
    nest::Node& t,
    nest::rport receptor_type,
    const CommonPropertiesType& )
  {
    ConnTestDummyNode dummy_target;
    ConnectionBase::check_connection_( dummy_target, s, t, receptor_type );
  }

  /**
   * Send an event to the receiver of this connection.
   * @param e The event to send
   * @param t Thread
   * @param t_lastspike Point in time of last spike sent.
   * @param cp Common properties to all synapses.
   */
  void update_dynamics();
  void send( nest::Event& e,
    nest::thread t,
    const CommonPropertiesType& cp );

  // The following methods contain mostly fixed code to forward the
  // corresponding tasks to corresponding methods in the base class and the w_
  // data member holding the weight.

  //! Store connection status information in dictionary
  void get_status( DictionaryDatum& d ) const;

  /**
   * Set connection status.
   *
   * @param d Dictionary with new parameter values
   * @param cm ConnectorModel is passed along to validate new delay values
   */
  void set_status( const DictionaryDatum& d, nest::ConnectorModel& cm );

  //! Allows efficient initialization on contstruction
  void
  set_weight( double w )
  {
    weight_ = w;
  }
};


template < typename targetidentifierT >
inline void
AibsConnection< targetidentifierT >::send( nest::Event& e,
  nest::thread t,

  const CommonPropertiesType& props )
{
  double del_t = e.get_stamp().get_ms() - t_lastspike_;
  nest::Node* target = get_target( t );

  double h=0.1; //currently, the time step for the Runge-Kutta 4th order is fixed at 0.1 msec. We can override it later. 
  double t0=0.0;

  int d = (int)(del_t / h);
  double n;
  double p;
  double Tau_r;
  double S;
  double p0;
    
    double kn1,  kn2,  kn3,  kn4;
    double kp1,  kp2,  kp3,  kp4;
    double ktr1, ktr2, ktr3, ktr4;
    double kS1,  kS2,  kS3,  kS4;
    double kp01, kp02, kp03, kp04;
    //std::cout<<e.get_stamp().get_ms()<<" "<<t_lastspike<<".\n";

   if (tau_r_<1e-10)
    {
        
        tau_r_=tau_r0_;
    }

    if (tau_f_<1e-10)
    {
        flag_f=0.0;
        p_=p0_; //p_ is fixed at the asymptotic value p0bar_
    }

    else
    {
        flag_f=1.0;
        p_=p0_;

    }

    if (tau_fdr_<1e-10)
    {
        flag_Tau_r=0.0;
	
    }

    else
    {
        flag_Tau_r=1.0;
        
    }

    if (tau_d_<1e-10)
    {
        flag_S=0.0;
    }

    else
    {
        flag_S=1.0;
    }

    if (tau_i_<1e-10)
    {
        flag_p0=0.0;
 
    }

    else
    {
        flag_p0=1.0;
        p0_=p0bar_;
        p_=p0bar_; //p_ is fixed at the asymptotic value p0bar_
    }
    //std::cout<<"f: "<<flag_f<<" Tau_r:"<<flag_Tau_r<<" S:"<<flag_S<<" p0:"<<flag_p0<<"\n";
    //std::cout<<"n: "<<n_<<" p:"<<p_<<" Tau_r:"<<tau_r_<<" S:"<<S_<<" p0:"<<p0_<<"\n";
    if (t_lastspike_>0.0001) //this ensures that the current spike is not the first spike. 
    {       
            //std::cout<<"update"<<"\n";
	    // initial impulse responses
            n = n_;
	    p=p_;
	    Tau_r=tau_r_;
	    S=S_;
	    p0=p0_;

            
	    for (int i=1; i<=d; i++)
	    {
		// Apply Runge Kutta Formulas to find
		// k1
		kn1   = h*dndt(t0, n, Tau_r);
		kp1   = h*dpdt(t0, p, p0, tau_f_)*flag_f;
		ktr1  = h*dtrdt(t0, Tau_r, tau_fdr_, tau_r0_)*flag_Tau_r;
		kS1   = h*dSdt(t0, S, tau_d_)*flag_S;
		kp01  = h*dp0dt(t0, p0, p0bar_,tau_i_)*flag_p0;
                //std::cout<<"kn1: "<<kn1<<" kp1:"<<kp1<<" Tau_r:"<<ktr1<<" kS1:"<<kS1<<" kp01:"<<kp01<<"\n";
	  
	 
		//k2
		kn2  = h*dndt(  t0 + 0.5*h,  n+0.5*kn1, Tau_r+0.5*ktr1);
		kp2  = h*dpdt(  t0 + 0.5*h,  p+0.5*kp1, p0+0.5*kp01, tau_f_)*flag_f;
		ktr2 = h*dtrdt( t0 + 0.5*h,  Tau_r + 0.5*ktr1, tau_fdr_, tau_r0_)*flag_Tau_r;
		kS2  = h*dSdt(  t0 + 0.5*h,  S+0.5*kS1, tau_d_)*flag_S;
		kp02 = h*dp0dt( t0 + 0.5*h,  p0+ 0.5*kp01, p0bar_,tau_i_)*flag_p0;

		//k3
		kn3  = h*dndt(  t0 + 0.5*h,  n+0.5*kn2, Tau_r+0.5*ktr2);
		kp3  = h*dpdt(  t0 + 0.5*h,  p+0.5*kp2, p0+0.5*kp02, tau_f_)*flag_f;
		ktr3 = h*dtrdt( t0 + 0.5*h,  Tau_r + 0.5*ktr2, tau_fdr_, tau_r0_)*flag_Tau_r;
		kS3  = h*dSdt(  t0 + 0.5*h,  S+0.5*kS2, tau_d_)*flag_S;
		kp03 = h*dp0dt( t0 + 0.5*h,  p0+ 0.5*kp02, p0bar_,tau_i_)*flag_p0;


		//k4
		kn4  = h*dndt(  t0 + h,  n+kn3, Tau_r+ktr3);
		kp4  = h*dpdt(  t0 + h,  p+kp3, p0+kp03, tau_f_)*flag_f;
		ktr4 = h*dtrdt( t0 + h,  Tau_r + ktr3, tau_fdr_, tau_r0_)*flag_Tau_r;
		kS4  = h*dSdt(  t0 + h,  S+kS3, tau_d_)*flag_S;
		kp04 = h*dp0dt( t0 + h,  p0+ kp03, p0bar_,tau_i_)*flag_p0;


	 
		// Update next value of y
		n = n + (1.0/6.0)*(kn1 + 2*kn2 + 2*kn3 + kn4);
		p = p + (1.0/6.0)*(kp1 + 2*kp2 + 2*kp3 + kp4);
		Tau_r = Tau_r + (1.0/6.0)*(ktr1 + 2*ktr2 + 2*ktr3 + ktr4);
		S = S + (1.0/6.0)*(kS1 + 2*kS2 + 2*kS3 + kS4);
		p0 = p0 + (1.0/6.0)*(kp01 + 2*kp02 + 2*kp03 + kp04);
		// Update next value of t
		t0 = t0 + h;
	 
	    }//end of for
          
	  n_=n;
	  p_=p;
	  S_=S;
	  p0_=p0;
	  tau_r_=Tau_r;
  }//end of if
  n_ = n_-p_*n_;
  p_=p_+p0_*(1-p_)*flag_f;
  tau_r_=tau_r_-a_fdr_*tau_r_*flag_Tau_r;
  S_=S_-a_D_*n_*p_*S_*flag_S;
  p0_=p0_-a_i_*p0_*flag_p0;

  // We use the current values for the spike number n.
  std::cout<<n_<<" "<<p_<<" "<<S_<<"\n";
  e.set_receiver( *target  );
  e.set_weight( p_ * n_ *S_*weight_ );
  // send the spike to the target
  e.set_delay( get_delay_steps() );
  e.set_rport( get_rport() );
  e();
  t_lastspike_ = e.get_stamp().get_ms();
}

template < typename targetidentifierT >
AibsConnection< targetidentifierT >::AibsConnection()
  : ConnectionBase()
  , weight_( 1.0 )
  , tau_f_(30.0)
  , tau_r_(2059.8)
  , tau_fdr_(569.8)
  , tau_r0_(2059.8)
  , tau_i_(100.0)
  , p0_(0.27)
  , n_(1.0)
  , p_(0.1)
  , a_fdr_(0.53)
  , a_D_(0.1)
  , a_i_(0.1)
  , S_(1.0)
  , p0bar_(0.5)  
  , tau_d_(30.0) 
{
}

template < typename targetidentifierT >
AibsConnection< targetidentifierT >::AibsConnection(
  const AibsConnection& rhs )
  : ConnectionBase( rhs )
  , weight_( rhs.weight_ )
  , tau_f_(rhs.tau_f_)
  , tau_r_(rhs.tau_r_)
  , tau_fdr_(rhs.tau_fdr_)
  , tau_r0_(rhs.tau_r0_)
  , tau_i_(rhs.tau_i_)
  , p0_(rhs.p0_)
  , n_(rhs.n_)
  , p_(rhs.p_)
  , a_fdr_(rhs.a_fdr_)
  , a_D_(rhs.a_D_)
  , a_i_(rhs.a_i_)
  , S_(rhs.S_)  
  , p0bar_(rhs.p0bar_)
  , tau_d_(rhs.tau_d_) 
{
}


template < typename targetidentifierT >
void
AibsConnection< targetidentifierT >::get_status(
  DictionaryDatum& d ) const
{
  ConnectionBase::get_status( d );

  def< double >( d, nest::names::weight, weight_ );


  def< double >( d, nest::names::tau_fac, tau_f_ );
  def< double >( d, nest::names::tau_rec, tau_r_ );
  def< double >( d, nest::names::tau_1, tau_fdr_ );
  def< double >( d, nest::names::tau_2, tau_r0_ );
  def< double >( d, nest::names::tau_eta, tau_i_ ); 
  def< double >( d, nest::names::x, p0_ );
  def< double >( d, nest::names::n, n_ );
  def< double >( d, nest::names::p, p_ );
  def< double >( d, nest::names::y_0, a_fdr_ );
  def< double >( d, nest::names::alpha_1, a_D_ );
  def< double >( d, nest::names::alpha_2, a_i_ );
  def< double >( d, nest::names::S, S_ );
  def< double >( d, nest::names::alpha, p0bar_);
  def< double >( d, nest::names::z, tau_d_ );
  def< long >( d, nest::names::size_of, sizeof( *this ) );
}

template < typename targetidentifierT >
void
AibsConnection< targetidentifierT >::set_status(
  const DictionaryDatum& d,
  nest::ConnectorModel& cm )
{
  ConnectionBase::set_status( d, cm );
  updateValue< double >( d, nest::names::weight, weight_ );

  updateValue< double >( d, nest::names::tau_fac, tau_f_ );
// if ( tau_fac_< 0.0 )
//  {
//    throw BadProperty( "Tau_fac must not be negative." );
//  }

  updateValue< double >( d, nest::names::tau_rec, tau_r_ );
//  if ( tau_rec_< 0.0 )
//  {
//    throw BadProperty( "Tau_rec must not be negative." );
//  }

  updateValue< double >( d, nest::names::tau_1, tau_fdr_ );
//  if ( tau_fdr_< 0.0 )
//  {
//    throw BadProperty( "Tau_fdr must not be negative." );
//  }

  updateValue< double >( d, nest::names::tau_2, tau_r0_ );
//  if ( tau_r0_< 0.0 )
//  {
//    throw BadProperty( "Tau_r0 must not be negative." );
//  }

  updateValue< double >( d, nest::names::tau_eta, tau_i_ );
//  if ( tau_i_< 0.0 )
//  {
//    throw BadProperty( "Tau_r0 must not be negative." );
//  }

  updateValue< double >( d, nest::names::x, p0_ );
//  if (p0_ > 1.0 || p0_ < 0.0 )
//  {
//    throw BadProperty( "p0 must be in [0,1]" );
// }  

  updateValue< double >( d, nest::names::n, n_ );
//  if (n_ > 1.0 || n_ < 0.0 )
//  {
//    throw BadProperty( "n must be in [0,1]" );
//  } 
  updateValue< double >( d, nest::names::p, p_ );
//  if (p_ > 1.0 || p_ < 0.0 )
//  {
//    throw BadProperty( "p must be in [0,1]" );
//  } 

  updateValue< double >( d, nest::names::y_0, a_fdr_ );
//  if (a_fdr_ > 1.0 || a_fdr_ < 0.0 )
//  {
//    throw BadProperty( "a_fdr must be in [0,1]" );
//  } 

  updateValue< double >( d, nest::names::alpha_1, a_D_ );
//  if (a_fdr_ > 1.0 || a_fdr_ < 0.0 )
//  {
//    throw BadProperty( "a_fdr must be in [0,1]" );
//  }

  updateValue< double >( d, nest::names::alpha_2, a_i_ );
//  if (a_fdr_ > 1.0 || a_fdr_ < 0.0 )
//  {
//    throw BadProperty( "a_fdr must be in [0,1]" );
//  } 


  updateValue< double >( d, nest::names::S, S_ );
//  if (S_ > 1.0 || S_ < 0.0 )
//  {
//    throw BadProperty( "S must be in [0,1]" );
//  } 

  updateValue< double >( d, nest::names::alpha, p0bar_ );
//  if (a_fdr_ > 1.0 || a_fdr_ < 0.0 )
//  {
//    throw BadProperty( "a_fdr must be in [0,1]" );
//  } 

  updateValue< double >( d, nest::names::z, tau_d_ );
//  if (tau_d_ > 1.0 || tau_d_ < 0.0 )
//  {
//    throw BadProperty( "tau_d must be in [0,1]" );
//  } 
}
} // namespace

#endif // drop_odd_spike_connection.h
