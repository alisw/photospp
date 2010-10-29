#ifndef _PhotosHepMCParticle_h_included_
#define _PhotosHepMCParticle_h_included_

/**
 * @class PhotosHepMCParticle
 *
 * @brief Interface to HepMC::GenParticle objects
 *
 * This class implements the virtual methods of
 * PhotosParticle. In this way it provides an
 * interface between the generic PhotosParticle class
 * and a HepMC::GenParticle object.
 *
 * @author Nadia Davidson
 * @date 17 June 2008
 */

#include <vector>

#include "HepMC/GenParticle.h"
#include "PhotosParticle.h"

using namespace std;

class PhotosHepMCParticle: public PhotosParticle{

 public:
  /** General constructor */
  PhotosHepMCParticle();

  /** Constructor which keeps a pointer to the HepMC::GenParticle*/
  PhotosHepMCParticle(HepMC::GenParticle * particle); 

  /** Constructor which creates a new HepMC::GenParticle and
       sets the properties pdg_id, statu and mass. */
  PhotosHepMCParticle(int pdg_id, int status, int mass);

  /** Destructor */
  ~PhotosHepMCParticle();
  
  /** return the HepMC::GenParticle */
  HepMC::GenParticle * getHepMC();

  /** Set the mothers of this particle via a vector of PhotosParticle*/
  void setMothers(std::vector<PhotosParticle*> mothers);

  /** Set the daughters of this particle via a vector of PhotosParticle*/
  void setDaughters(std::vector<PhotosParticle*> daughters);

  /** Add a new daughter to the end vertex of this particle */ 
  void addDaughter(PhotosParticle* daughter);

  /** Returns the mothers of this particle via a vector of PhotosParticle */
  std::vector<PhotosParticle*> getMothers();

  /** Returns the daughters of this particle via a vector of PhotosParticle */
  std::vector<PhotosParticle*> getDaughters();

  /** Set the PDG ID code of this particle */
  void setPdgID(int pdg_id);

  /** Set the status of this particle */
  void setStatus(int statu);

  /** Set the mass of this particle */
  void setMass(double mass);

  /** Get the PDG ID code of this particle */
  int getPdgID();

  /** Get the status of this particle */
  int getStatus();

  /** Get the barcode of this particle */
  int getBarcode();

  /** check that the 4 momentum in conserved at the vertices producing
      and ending this particle */
  bool checkMomentumConservation();

  /** Create a new particle of type PhotosHepMCParticle, with the given
      properties. The new particle bares no relations to this
      particle, but it provides a way of creating a instance of
      this derived class. eg. createNewParticle() is used inside
      filhep_() so that a PhotosHepMCParticle can be created without
      the method having explicit knowledge of the PhotosHepMCParticle 
      class */
  PhotosHepMCParticle * createNewParticle(int pdg_id, int status, double mass,
				       double px, double py,
				       double pz, double e);

  /** Print some information about this particle to standard output */
  void print();

  /** Returns the px component of the four vector*/
  double getPx();

  /** Returns the py component of the four vector */
  double getPy();

  /** Returns the pz component of the four vector */
  double getPz();

  /** Returns the energy component of the four vector */
  double getE();

  /** Returns the mass taken from event record */
  double getMass();

  /** Set the px component of the four vector */
  void setPx( double px );

  /** Set the px component of the four vector */
  void setPy( double py );

  /** Set the pz component of the four vector */
  void setPz( double pz );

  /** Set the energy component of the four vector */
  void setE( double e );

 private:
  /** Internal function used to clear particles from the vector */
  void clear(std::vector<PhotosParticle*> v);

  /** A pointer to the HepMC::GenParticle particle */
  HepMC::GenParticle * m_particle;

  /** A vector of this particles mothers */
  std::vector<PhotosParticle*> m_mothers;

  /** A vector of this particles daughters */
  std::vector<PhotosParticle*> m_daughters;

  /** list to keep track of new particles which have been
      created from this one, so we can call their destructor later */
  std::vector<PhotosParticle*> m_created_particles;

};

#endif  

