/*! \file ECModels.h
    \brief Contains declaration of classes that implements models for fields from eddy currents.

    \author Jesper Andersson
    \version 1.0b, Sep., 2012.
*/
// Declarations of classes that implements a hirearchy
// of models for fields from eddy currents induced by
// diffusion gradients.
// 
// ECModels.h
//
// Jesper Andersson, FMRIB Image Analysis Group
//
// Copyright (C) 2011 University of Oxford 
//

#ifndef ECModels_h
#define ECModels_h

#include <cstdlib>
#include <string>
#include <vector>
#include <cmath>
#include <boost/shared_ptr.hpp>
#include "newmat.h"
#include "newimage/newimageall.h"
#include "miscmaths/miscmaths.h"
#include "topup/topup_file_io.h"
#include "EddyHelperClasses.h"

namespace EDDY {

/****************************************************************//**
*
* \brief Virtual base class for classes used to model the fields
* that may result from eddy currents.
*
* The classes in this hierarchy manages eddy current (EC) parameters
* and movement parameters for one scan. We can set the parameters 
* with one call and obtain the resulting field with another. By
* deriving a set of classes from a virtual base class we are able
* to use the same code to estimate the parameters for different 
* EC-models.
*
********************************************************************/ 
class ScanECModel
{
public:
  ScanECModel() {}
  ScanECModel(const NEWMAT::ColumnVector& mp,
              const NEWMAT::ColumnVector& ep) : _mp(mp), _ep(ep) {}
  virtual ~ScanECModel() {}

  /// Indicates if a field offset is modeled or not.
  virtual bool HasFieldOffset() const = 0;
  /// Returns the field offset.
  virtual double GetFieldOffset() const = 0;
  /// Set field offset.
  virtual void SetFieldOffset(double ofst) = 0;

  // Returns linear parameters used for modelling a field offset
  virtual NEWMAT::RowVector GetLinearParameters() const { NEWMAT::RowVector empty; return(empty); }

  // Return the total number of parameters
  unsigned int NParam(Parameters wp=ALL) const { 
    switch(wp) { case MOVEMENT: return(_mp.Nrows()); case EC: return(_ep.Nrows()); case ALL: return(_mp.Nrows() + _ep.Nrows()); }
    return(_mp.Nrows() + _ep.Nrows()); // To silence warning
  }

  /// Get all parameters for the category given by mp.
  NEWMAT::ColumnVector GetParams(Parameters wp=ALL) const { 
    switch(wp) { case MOVEMENT: return(_mp); case EC: return(_ep); case ALL: return(_mp & _ep); }
    return(_mp & _ep);  // To silence warning
  }

  /// Get parameter given by indx, where indx numbers parameters starting with the movement parameters
  double GetParam(unsigned int indx, Parameters wp=ALL) const {
    try { switch(wp) {
      case MOVEMENT: return(mov_index2ref(indx));
      case EC: return(ec_index2ref(indx)); 
      case ALL: return(index2ref(indx)); } }
    catch(...) { throw EddyException("ScanECModel::GetParam: Caught exception"); }
    return(0.0); // To silence warning
  }    

  // Set all parameters in specified category
  void SetParams(const NEWMAT::ColumnVector& p, Parameters wp=ALL) {
    switch (wp) {
    case MOVEMENT: 
      if (p.Nrows() != _mp.Nrows()) throw EddyException("ScanECModel::SetParams: Wrong number of parameters"); else _mp = p; break;
    case EC:
      if (p.Nrows() != _ep.Nrows()) throw EddyException("ScanECModel::SetParams: Wrong number of parameters"); else _ep = p; break;
    case ALL:
      if (p.Nrows() != (_mp.Nrows()+_ep.Nrows())) throw EddyException("ScanECModel::SetParams: Wrong number of parameters");
      else _mp = p.Rows(1,_mp.Nrows()); _ep = p.Rows(_mp.Nrows()+1,p.Nrows());
      break;
    }
  }

  // Set parameter given by indx, where indx numbers all parameters
  // zero-offset starting with the movement parameters
  void SetParam(unsigned int indx, double p, Parameters wp=ALL) {
    try { switch(wp) {
      case MOVEMENT: mov_index2ref(indx) = p; break; 
      case EC: ec_index2ref(indx) = p; break; 
      case ALL: index2ref(indx) = p; break; } }
    catch(...) { throw EddyException("ScanECModel::SetParam: Caught exception"); }
  }
  // Returns matrix denoted \mathbf{R} in paper.
  NEWMAT::Matrix ForwardMovementMatrix(const NEWIMAGE::volume<float>& scan) const { return(TOPUP::MovePar2Matrix(_mp,scan)); }
  // Returns matrix denoted \nathbf{R}^{-1} in paper.
  NEWMAT::Matrix InverseMovementMatrix(const NEWIMAGE::volume<float>& scan) const { return(ForwardMovementMatrix(scan).i()); }

  // Return the number of parameters that are updated as part of the estimation.
  // For some models this may be different from total number of parameters
  virtual unsigned int NDerivs(Parameters wp=ALL) const { return(NParam(wp)); }

  // The following get/set routines indexes parameters from
  // 0 - NDerivs()-1  i.e. it ignores any parameters that are
  // not estimated for the particular model.
  virtual double GetDerivParam(unsigned int dindx, Parameters wp=ALL) const {  
    try { return(GetParam(dindx,wp)); }
    catch(...) { throw EddyException("ScanECModel::GetDerivParam: Caught exception"); }
  }

  virtual void SetDerivParam(unsigned int dindx, double p, Parameters wp=ALL) {
    try { SetParam(dindx,p,wp); }
    catch(...) { throw EddyException("ScanECModel::SetDerivParam: Caught exception"); }
  }

  // Returns a pointer to a clone of itself
  virtual boost::shared_ptr<ScanECModel> Clone() const = 0;

  // Returns sutiable scales for evaluating numerical derivatives
  virtual double GetDerivScale(unsigned int dindx, Parameters wp=ALL) const = 0;

  // Return eddy current-induced field in Hz. Denoted e(\mathbf{b}) in paper.
  virtual NEWIMAGE::volume<float> ECField(const NEWIMAGE::volume<float>& scan) const = 0;

protected:
  NEWMAT::ColumnVector _mp; // dx dy dz rx ry rz, N.B. different from MJ
  NEWMAT::ColumnVector _ep; // Different for different models

  double index2ref(unsigned int indx) const { 
    if (indx < static_cast<unsigned int>(_mp.Nrows())) return(_mp(indx+1));
    else if (indx < static_cast<unsigned int>(_mp.Nrows()+_ep.Nrows())) return(_ep(indx-_mp.Nrows()+1));
    else throw EddyException("ScanECModel::index2ref const: Index out of range");
  }
  double& index2ref(unsigned int indx) { 
    if (indx < static_cast<unsigned int>(_mp.Nrows())) return(_mp(indx+1));
    else if (indx < static_cast<unsigned int>(_mp.Nrows()+_ep.Nrows())) return(_ep(indx-_mp.Nrows()+1));
    else throw EddyException("ScanECModel::index2ref: Index out of range");
  }
  double mov_index2ref(unsigned int indx) const { 
    if (indx < static_cast<unsigned int>(_mp.Nrows())) return(_mp(indx+1));
    else throw EddyException("ScanECModel::mov_index2ref const: Index out of range");
  }
  double& mov_index2ref(unsigned int indx) { 
    if (indx < static_cast<unsigned int>(_mp.Nrows())) return(_mp(indx+1));
    else throw EddyException("ScanECModel::mov_index2ref: Index out of range");
  }
  double ec_index2ref(unsigned int indx) const { 
    if (indx < static_cast<unsigned int>(_ep.Nrows())) return(_ep(indx+1));
    else throw EddyException("ScanECModel::ec_index2ref const: Index out of range");
  }
  double& ec_index2ref(unsigned int indx) { 
    if (indx < static_cast<unsigned int>(_ep.Nrows())) return(_ep(indx+1));
    else throw EddyException("ScanECModel::ec_index2ref: Index out of range");
  }
};

//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//
// Class PolynomialScanECModel (Polynomial Scan Eddy Current Model)
//
// A virtual base class for polynomial models (linear, quadratic etc).
// 
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

class PolynomialScanECModel : public ScanECModel
{
public:
  PolynomialScanECModel(bool field=false) : ScanECModel() {
    _mp.ReSize(6); _mp=0.0; _nepd = 0;
  }
  PolynomialScanECModel(const NEWMAT::ColumnVector& mp, const NEWMAT::ColumnVector& ep, bool field=false) : ScanECModel(mp,ep), _nepd(0) {}
  virtual ~PolynomialScanECModel() {} 
  virtual bool HasFieldOffset() const = 0;
  virtual double GetFieldOffset() const { if (HasFieldOffset()) return(_ep(_nepd)); else return(0.0); }
  virtual void SetFieldOffset(double ofst) { if (HasFieldOffset()) _ep(_nepd) = ofst; }
  virtual NEWMAT::RowVector GetLinearParameters() const { return(_ep.Rows(1,3).t()); }

  // Return the number of parameters that are updated as part of the estimation.
  // For some models this may be different from total number of parameters
  virtual unsigned int NDerivs(Parameters wp=ALL) const { 
    switch(wp) { case MOVEMENT: return(_mp.Nrows()); case EC: return(_nepd); case ALL: return(_mp.Nrows()+_nepd); } 
    return(_mp.Nrows()+_nepd); // To silence warning
  }
  // The following get/set routines indexes parameters from
  // 0 - NDerivs()-1  i.e. it ignores any parameters that are
  // not estimated for the particular model.
  virtual double GetDerivParam(unsigned int dindx, Parameters wp=ALL) const {  
    try { switch(wp) { case MOVEMENT: return(mov_dindex2ref(dindx)); case EC: return(ec_dindex2ref(dindx)); case ALL: return(dindex2ref(dindx)); } }
    catch(...) { throw EddyException("PolynomialScanECModel::GetDerivParam: Caught exception"); }
    return(0.0); // To silence warning
  }
  virtual void SetDerivParam(unsigned int dindx, double p, Parameters wp=ALL) {
    try { switch(wp) { 
      case MOVEMENT: mov_dindex2ref(dindx) = p; break;
      case EC: ec_dindex2ref(dindx) = p; break; 
      case ALL: dindex2ref(dindx) = p; } }
    catch(...) { throw EddyException("PolynomialScanECModel::SetDerivParam: Caught exception"); }
  }
  // Returns suitable scales for evaluating numerical derivatives
  virtual double GetDerivScale(unsigned int dindx, Parameters wp=ALL) const { 
    try { switch(wp) {
      case MOVEMENT: return(get_mov_deriv_scale(dindx));
      case EC: return(get_ec_deriv_scale(dindx-_mp.Nrows()));
      case ALL: if (dindx < static_cast<unsigned int>(_mp.Nrows())) return(get_mov_deriv_scale(dindx)); else return(get_ec_deriv_scale(dindx-_mp.Nrows())); } }
    catch(...) { throw EddyException("PolynomialScanECModel::GetDerivScale: Caught Exception"); }
    return(0.0); // To silence warning
  }

  virtual boost::shared_ptr<ScanECModel> Clone() const = 0; 

  // Return eddy current-induced field in Hz. Denoted e(\mathbf{b}) in paper.
  virtual NEWIMAGE::volume<float> ECField(const NEWIMAGE::volume<float>& scan) const = 0;

protected:
  unsigned int _nepd;  // Number of Eddy Parameter Derivatives (might depend on if field was set)

private:

  double dindex2ref(unsigned int indx) const {
    if (indx < static_cast<unsigned int>(_mp.Nrows())) return(_mp(indx+1));
    else if (indx < static_cast<unsigned int>(_mp.Nrows())+_nepd) return(_ep(indx-_mp.Nrows()+1));
    else throw EddyException("PolynomialScanECModel::dindex2ref const: Index out of range");
  }
  double& dindex2ref(unsigned int indx) {
    if (indx < static_cast<unsigned int>(_mp.Nrows())) return(_mp(indx+1));
    else if (indx < static_cast<unsigned int>(_mp.Nrows())+_nepd) return(_ep(indx-_mp.Nrows()+1));
    else throw EddyException("PolynomialScanECModel::dindex2ref: Index out of range");
  }
  double mov_dindex2ref(unsigned int indx) const {
    if (indx < static_cast<unsigned int>(_mp.Nrows())) return(_mp(indx+1));
    else throw EddyException("PolynomialScanECModel::mov_dindex2ref const: Index out of range");
  }
  double& mov_dindex2ref(unsigned int indx) {
    if (indx < static_cast<unsigned int>(_mp.Nrows())) return(_mp(indx+1));
    else throw EddyException("PolynomialScanECModel::mov_dindex2ref: Index out of range");
  }
  double ec_dindex2ref(unsigned int indx) const { 
    if (indx < _nepd) return(_ep(indx+1));
    else throw EddyException("PolynomialScanECModel::ec_dindex2ref const: Index out of range");
  }
  double& ec_dindex2ref(unsigned int indx) { 
    if (indx < _nepd) return(_ep(indx+1));
    else throw EddyException("PolynomialScanECModel::ec_dindex2ref: Index out of range");
  }
  virtual double get_mov_deriv_scale(unsigned int dindx) const = 0;
  virtual double get_ec_deriv_scale(unsigned int dindx) const = 0;
};

//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//
// Class LinearScanECModel (Linear Scan Eddy Current Model)
//
// This class models the eddy curents as a linear combination of
// linear gradients in the x-, y- and z-directions. The assumption
// behind this is that the eddy currents resides mainly in the 
// gradient coils and also that the gradient coils are close to
// linear.
// The _ep field contains:
// dfdx (shear), dfdy (zoom), dfdz (z-dependent trans) all in Hz/mm and df (trans in Hz)
// Note that the translations (e.g. dfdx (EC gradient in x-direction) to shear) 
// only makes sense if the phase-encode is in the y-direction.
// The df field is there to model any difference between the centre
// of the FOV and the iso-centre of the scanner.
//
// If no susceptibility off-resonance field is specified df is strictly
// speaking redundant since it is identical to a subject movement
// (y-translation). It is however useful to retain as a parameter
// for the sake of modelling EC parameters as a function of 
// diffusion gradients. In the updates it will not have a derivative
// and will hence not be updated. It will instead be set by the 
// higher level modelling of the parameters.
//
// If a susceptibility off-resonance field is specified df will enter
// into the transform compared to subject y-translation and it should
// be possible to directly separate the two. In this case it will
// have a derivative and will be updated as part of the first level
// estimation.
//
// {{{ @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

class LinearScanECModel : public PolynomialScanECModel
{
public:
  LinearScanECModel(bool field=false) : PolynomialScanECModel() {
    _mp.ReSize(6); _mp=0.0; _ep.ReSize(4); _ep=0.0;
    _nepd = (field) ? 4 : 3;
  }
  LinearScanECModel(const NEWMAT::ColumnVector& mp, const NEWMAT::ColumnVector& ep, bool field=false) : PolynomialScanECModel(mp,ep) {
    if (_mp.Nrows() != 6) throw EddyException("LinearScanECModel: mp must have 6 elements");
    _nepd = (field) ? 4 : 3;
    if (_ep.Nrows() != int(_nepd)) throw EddyException("LinearScanECModel: Wrong number of elements for ep");
  }
  virtual ~LinearScanECModel() {} 
  bool HasFieldOffset() const { return(_nepd==4); }

  virtual boost::shared_ptr<ScanECModel> Clone() const { return(boost::shared_ptr<ScanECModel>( new LinearScanECModel(*this))); }

  // Return eddy current-induced field in Hz. Denoted e(\mathbf{b}) in paper.
  virtual NEWIMAGE::volume<float> ECField(const NEWIMAGE::volume<float>& scan) const;

private:

  double get_mov_deriv_scale(unsigned int dindx) const { 
    if (dindx < 3) return(1e-3); 
    else if (dindx < 6) return(1e-5); 
    else throw EddyException("LinearScanECModel::get_mov_deriv_scale: Index out of range");
  }
  double get_ec_deriv_scale(unsigned int dindx) const { 
    if (dindx < 3) return(1e-3);
    else if (dindx < _nepd) return(1e-2); 
    else throw EddyException("LinearScanECModel::get_ec_deriv_scale: Index out of range");
  }
};

//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//
// Class QuadraticScanECModel (Quadratic Scan Eddy Current Model)
//
// This class models the eddy curents as a second order polynomial
// expansion of  combination of gradients in the x-, y- and z-directions. 
//
// The _ep field contains:
// dfdx (shear), dfdy (zoom), dfdz (z-dependent trans) all in Hz/mm 
// followed by dfdx^2, dfdy^2, dfdz^2, dfdx*dfdy, dfdx*dfdz and dfdy*dfdz 
// and finally df (trans in Hz).
// The quadratic components are (arbitrarily) scaled to have ~the same value (in Hz) 
// at the edge of the FOV as has the linear terms. This is done to ensure an
// update matrix with reasonable condition number.
// Note that the translations (e.g. dfdx (EC gradient in x-direction) to shear) 
// only makes sense if the phase-encode is in the y-direction.
// The df field is there to model any difference between the centre
// of the FOV and the iso-centre of the scanner.
//
// If no susceptibility off-resonance field is specified df is strictly
// speaking redundant since it is identical to a subject movement
// (y-translation). It is however useful to retain as a parameter
// for the sake of modelling EC parameters as a function of 
// diffusion gradients. In the updates it will not have a derivative
// and will hence not be updated. It will instead be set by the 
// higher level modelling of the parameters.
//
// If a susceptibility off-resonance field is specified df will enter
// into the transform compared to subject y-translation and it should
// be possible to directly separate the two. In this case it will
// have a derivative and will be updated as part of the first level
// estimation.
//
// {{{ @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

class QuadraticScanECModel : public PolynomialScanECModel
{
public:
  QuadraticScanECModel(bool field=false) : PolynomialScanECModel() {
    _mp.ReSize(6); _mp=0.0; _ep.ReSize(10); _ep=0.0;
    _nepd = (field) ? 10 : 9;
  }
  QuadraticScanECModel(const NEWMAT::ColumnVector& mp, const NEWMAT::ColumnVector& ep, bool field=false) : PolynomialScanECModel(mp,ep) {
    if (_mp.Nrows() != 6) throw EddyException("QuadraticScanECModel: mp must have 6 elements");
    _nepd = (field) ? 10 : 9;
    if (_ep.Nrows() != int(_nepd)) throw EddyException("QuadraticScanECModel: Wrong number of elements for ep");
  }
  virtual ~QuadraticScanECModel() {} 
  bool HasFieldOffset() const { return(_nepd==10); }

  boost::shared_ptr<ScanECModel> Clone() const { return(boost::shared_ptr<ScanECModel>( new QuadraticScanECModel(*this))); }

  // Return eddy current-induced field in Hz. Denoted e(\mathbf{b}) in paper.
  NEWIMAGE::volume<float> ECField(const NEWIMAGE::volume<float>& scan) const;

private:
  double get_mov_deriv_scale(unsigned int dindx) const { 
    if (dindx < 3) return(1e-3); 
    else if (dindx < 6) return(1e-5); 
    else throw EddyException("QuadraticScanECModel::get_mov_deriv_scale: Index out of range");
  }
  double get_ec_deriv_scale(unsigned int dindx) const { 
    if (dindx < 3) return(1e-3);
    else if (dindx < 9) return(1e-5);
    else if (dindx < _nepd) return(1e-2); 
    else throw EddyException("QuadraticScanECModel::get_ec_deriv_scale: Index out of range");
  }
};

//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//
// Class CubicScanECModel (Cubic Scan Eddy Current Model)
//
// This class models the eddy curents as a third order polynomial
// expansion of  combination of gradients in the x-, y- and z-directions. 
//
// The _ep field contains:
// dfdx (shear), dfdy (zoom), dfdz (z-dependent trans) all in Hz/mm 
// followed by dfdx^2, dfdy^2, dfdz^2, dfdx*dfdy, dfdx*dfdz and dfdy*dfdz, 
// followed by dfdx^3, dfdy^3, dfdz^3, dfdx^2*dfdy, dfdx^2*dfdz,
// dfdy^2*dfdx, dfdy^2*dfdz, dfdz^2*dfdx, dfdz^2*dfdy, dfdx*dfdy*dfdz
//  and finally df (trans in Hz).
// The quadratic and cubic components are (arbitrarily) scaled to have ~the 
// same value (in Hz) at the edge of the FOV as has the linear terms. 
// This is done to ensure an update matrix with reasonable condition number.
// Note that the translations (e.g. dfdx (EC gradient in x-direction) to shear) 
// only makes sense if the phase-encode is in the y-direction.
// The df field is there to model any difference between the centre
// of the FOV and the iso-centre of the scanner.
//
// If no susceptibility off-resonance field is specified df is strictly
// speaking redundant since it is identical to a subject movement
// (y-translation). It is however useful to retain as a parameter
// for the sake of modelling EC parameters as a function of 
// diffusion gradients. In the updates it will not have a derivative
// and will hence not be updated. It will instead be set by the 
// higher level modelling of the parameters.
//
// If a susceptibility off-resonance field is specified df will enter
// into the transform compared to subject y-translation and it should
// be possible to directly separate the two. In this case it will
// have a derivative and will be updated as part of the first level
// estimation.
//
// {{{ @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

class CubicScanECModel : public PolynomialScanECModel
{
public:
  CubicScanECModel(bool field=false) : PolynomialScanECModel() {
    _mp.ReSize(6); _mp=0.0; _ep.ReSize(20); _ep=0.0;
    _nepd = (field) ? 20 : 19;
  }
  CubicScanECModel(const NEWMAT::ColumnVector& mp, const NEWMAT::ColumnVector& ep, bool field=false) : PolynomialScanECModel(mp,ep) {
    if (_mp.Nrows() != 6) throw EddyException("CubicScanECModel: mp must have 6 elements");
    _nepd = (field) ? 20 : 19;
    if (_ep.Nrows() != int(_nepd)) throw EddyException("CubicScanECModel: Wrong number of elements for ep");
  }
  virtual ~CubicScanECModel() {} 
  bool HasFieldOffset() const { return(_nepd==20); }

  boost::shared_ptr<ScanECModel> Clone() const { return(boost::shared_ptr<ScanECModel>( new CubicScanECModel(*this))); }

  // Return eddy current-induced field in Hz. Denoted e(\mathbf{b}) in paper.
  NEWIMAGE::volume<float> ECField(const NEWIMAGE::volume<float>& scan) const;

private:
  double get_mov_deriv_scale(unsigned int dindx) const { 
    if (dindx < 3) return(1e-3); 
    else if (dindx < 6) return(1e-5); 
    else throw EddyException("CubicScanECModel::get_mov_deriv_scale: Index out of range");
  }
  double get_ec_deriv_scale(unsigned int dindx) const { 
    if (dindx < 3) return(1e-3);
    else if (dindx < 9) return(1e-5);
    else if (dindx < 19) return(1e-7);
    else if (dindx < _nepd) return(1e-2); 
    else throw EddyException("CubicScanECModel::get_ec_deriv_scale: Index out of range");
  }
};

//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//
// Class MovementScanECModel (Movement Scan Eddy Current Model)
//
// This class doesn't model the eddy curents at all, and simply 
// uses the rigid body model. This is done to create a polymorphism
// so we can use the same basic code for the b0 scans as for the
// diffusion weighted ones.
//
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

class MovementScanECModel : public ScanECModel
{
public:
  MovementScanECModel() : ScanECModel() {
    _mp.ReSize(6); _mp=0.0;
  }
  MovementScanECModel(const NEWMAT::ColumnVector& mp) {
    if (_mp.Nrows() != 6) throw EddyException("MovementScanECModel: mp must have 6 elements");
    _mp = mp;
  }
  virtual ~MovementScanECModel() {} 
  virtual bool HasFieldOffset() const { return(false); }
  virtual double GetFieldOffset() const { return(0.0); }
  virtual void SetFieldOffset(double ofst) { }

  // Return the number of parameters that are updated as part of the estimation.
  // For some models this may be different from total number of parameters
  unsigned int NDerivs(Parameters wp=ALL) const { 
    if (wp == EC) throw EddyException("MovementScanECModel::NDerivs: Model has no EC parameters");
    else return(_mp.Nrows());
  }
  // The following get/set routines indexes parameters from
  // 0 - NDerivs()-1  i.e. it ignores any parameters that are
  // not estimated for the particular model.
  double GetDerivParam(unsigned int dindx, Parameters wp=ALL) const {
    if (wp == EC) throw EddyException("MovementScanECModel::GetDerivParam: Model has no EC parameters");
    if (dindx > 5) throw EddyException("MovementScanECModel::GetDerivParam: index out of range");  
    return(_mp(dindx+1));
  }
  void SetDerivParam(unsigned int dindx, double p, Parameters wp=ALL) {
    if (wp == EC) throw EddyException("MovementScanECModel::SetDerivParam: Model has no EC parameters");
    if (dindx > 5) throw EddyException("MovementScanECModel::SetDerivParam: index out of range");  
    _mp(dindx+1) = p; // To silence warning
  }
  // Returns suitable scales for evaluating numerical derivatives
  double GetDerivScale(unsigned int dindx, Parameters wp=ALL) const { 
    if (wp == EC) throw EddyException("MovementScanECModel::GetDerivScale: Model has no EC parameters");
    if (dindx > 5) throw EddyException("MovementScanECModel::GetDerivScale: index out of range");  
    return(get_mov_deriv_scale(dindx));
  }

  boost::shared_ptr<ScanECModel> Clone() const { return(boost::shared_ptr<ScanECModel>( new MovementScanECModel(*this))); }

  // Return eddy current-induced field in Hz. Denoted e(\mathbf{b}) in paper.
  NEWIMAGE::volume<float> ECField(const NEWIMAGE::volume<float>& scan) const;

private:
  double get_mov_deriv_scale(unsigned int dindx) const { 
    if (dindx < 3) return(1e-3); 
    else if (dindx < 6) return(1e-5); 
    else throw EddyException("MovementScanECModel::get_mov_deriv_scale: Index out of range");
  }
};

} // End namespace EDDY

#endif // End #ifndef ECModels_h
