#include "PhotosHEPEVTParticle.h"
#include "Log.h"
typedef Photos::Log Log;

PhotosHEPEVTParticle::~PhotosHEPEVTParticle()
{
  // Cleanup particles that do not belong to event
  for(unsigned int i=0;i<cache.size();i++)
    if(cache[i]->m_barcode<0)
      delete cache[i];
}

PhotosHEPEVTParticle::PhotosHEPEVTParticle(int pdgid, int status, double px, double py, double pz, double e, double m, int ms, int me, int ds, int de){
  m_px = px;
  m_py = py;
  m_pz = pz;
  m_e  = e;
  m_generated_mass = m;

  m_pdgid  = pdgid;
  m_status = status;

  m_mother_start   = ms;
  m_mother_end     = me;
  m_daughter_start = ds;
  m_daughter_end   = de;

  m_barcode = -1;
  m_event = NULL;
  
  isInRestFrame = false;
}

/** Add a new daughter to this particle */
void PhotosHEPEVTParticle::addDaughter(PhotosParticle* daughter)
{
  Log::Info()<<"Ha, Haaaaa... not happening!"<<endl;
}

void PhotosHEPEVTParticle::setMothers(vector<PhotosParticle*> mothers){

  int beg = 65535, end = -1;

  // If this particle has not yet been added to the event record
  // then add it to the mothers' event record
  if(m_barcode<0 && mothers.size()>0)
  {
    PhotosHEPEVTEvent *evt = ((PhotosHEPEVTParticle*)mothers[0])->m_event;
    evt->addParticle(this);
  }

  for(unsigned int i=0;i<mothers.size();i++)
  {
    int bc = mothers[i]->getBarcode();
    if(bc<0) Log::Fatal("Fatal..\n");
    if(bc<beg) beg = bc;
    if(bc>end) end = bc;
  }
  if(end == -1) beg = -1;

  m_mother_start = beg;
  m_mother_end   = end;

}

void PhotosHEPEVTParticle::setDaughters(vector<PhotosParticle*> daughters){

  int beg = 65535, end = -1;

  for(unsigned int i=0;i<daughters.size();i++)
  {
    int bc = daughters[i]->getBarcode();
    if(bc<0) Log::Fatal("PhotosHEPEVTParticle::setDaughters: all daughters has to be in event record first");
    if(bc<beg) beg = bc;
    if(bc>end) end = bc;
  }
  if(end == -1) beg = -1;

  m_daughter_start = beg;
  m_daughter_end   = end;
}

std::vector<PhotosParticle*> PhotosHEPEVTParticle::getMothers(){

  std::vector<PhotosParticle*> mothers;

  if(m_event && m_mother_end>=0)
  {
    for(int i=m_mother_start;i<=m_mother_end;i++)
    {
      PhotosParticle *p = m_event->getParticle(i);
      if( p==NULL )
      {
        Log::Warning()<<"PhotosHEPEVTParticle::getMothers(): No particle with index "<<i<<endl;
        continue;
      }

      mothers.push_back(p);
    }
  }
  return mothers;
}

std::vector<PhotosParticle*> PhotosHEPEVTParticle::getDaughters(){

  std::vector<PhotosParticle*> daughters;

  if(!m_event) return daughters;

  // Check if m_daughter_start and m_daughter_end are set
  // If not - try to set get daughters list from event
  if(m_daughter_end<0)
  {
    int min_d=65535, max_d=-1;
    for(int i=0;i<m_event->getParticleCount();i++)
    {
      if(m_event->getParticle(i)->isDaughterOf(this))
      {
        if(i<min_d) min_d = i;
        if(i>max_d) max_d = i;
      }
    }
    if(max_d>=0)
    {
      m_daughter_start = min_d;
      m_daughter_end   = max_d;
      m_status         = 2;
    }
  }

  if(m_daughter_end>=0)
  {
    for(int i=m_daughter_start;i<=m_daughter_end;i++)
    {
      PhotosParticle *p = m_event->getParticle(i);
      if(p==NULL)
      {
        Log::Warning()<<"PhotosHEPEVTParticle::getDaughters(): No particle with index "<<i<<endl;
        return daughters;
      }

      daughters.push_back(p);
    }
  }

  return daughters;
}

bool PhotosHEPEVTParticle::checkMomentumConservation(){

  if(!m_event)           return true;
  if(m_daughter_end < 0) return true;

  PhotosHEPEVTParticle *buf = m_event->getParticle(m_daughter_start);

  int mother_start = buf->m_mother_start;
  int mother_end   = buf->m_mother_end;

  double px =0.0, py =0.0, pz =0.0;

  for(int i=m_daughter_start;i<=m_daughter_end;i++)
  {
    buf = m_event->getParticle(i);
    px += buf->getPx();
    py += buf->getPy();
    pz += buf->getPz();
  }

  for(int i=mother_start;i<=mother_end;i++)
  {
    buf = m_event->getParticle(i);
    px -= buf->getPx();
    py -= buf->getPy();
    pz -= buf->getPz();
  }

  if( sqrt( px*px + py*py + pz*pz ) > Photos::momentum_conservation_threshold )
  {
    Log::RedirectOutput( Log::Warning()<<"Momentum not conserved in vertex: " );
    for(int i=mother_start;    i<=mother_end;    i++) m_event->getParticle(i)->print();
    for(int i=m_daughter_start;i<=m_daughter_end;i++) m_event->getParticle(i)->print();
    Log::RevertOutput();
    return false;
  }

  return true;
}

PhotosHEPEVTParticle * PhotosHEPEVTParticle::createNewParticle(
                        int pdg_id, int status, double mass,
                        double px, double py, double pz, double e){
  cache.push_back(new PhotosHEPEVTParticle(pdg_id,status,px,py,pz,e,mass,-1,-1,-1,-1));
  return cache.back();
}

bool PhotosHEPEVTParticle::isDaughterOf(PhotosHEPEVTParticle *p)
{
  int bc = p->getBarcode();
  if(bc>=m_mother_start && bc<=m_mother_end) return true;

  return false;
}

bool PhotosHEPEVTParticle::isMotherOf  (PhotosHEPEVTParticle *p)
{
  int bc = p->getBarcode();
  if(bc>=m_daughter_start && bc<=m_daughter_end) return true;

  return false;
}

void PhotosHEPEVTParticle::print(){
  printf("P: (%2i) %6i %2i | %11.4e %11.4e %11.4e %11.4e | %11.4e | M: %2i %2i | D: %2i %2i\n",
          m_barcode, m_pdgid, m_status, m_px, m_py, m_pz, m_e, m_generated_mass,
          m_mother_start, m_mother_end,   m_daughter_start, m_daughter_end);
}

/******** Getter and Setter methods: ***********************/

void PhotosHEPEVTParticle::setPdgID(int pdg_id){
  m_pdgid = pdg_id;
}

void PhotosHEPEVTParticle::setStatus(int status){
  m_status = status;
}

void PhotosHEPEVTParticle::setMass(double mass){
  m_generated_mass = mass;
}

int PhotosHEPEVTParticle::getPdgID(){
  return m_pdgid;
}

int PhotosHEPEVTParticle::getStatus(){
  return m_status;
}

double PhotosHEPEVTParticle::getMass(){
  return m_generated_mass;
}

inline double PhotosHEPEVTParticle::getPx(){
  return m_px;
}

inline double PhotosHEPEVTParticle::getPy(){
  return m_py;
}

double PhotosHEPEVTParticle::getPz(){
  return m_pz;
}

double PhotosHEPEVTParticle::getE(){
  return m_e;
}

void PhotosHEPEVTParticle::setPx(double px){
  m_px = px;
}

void PhotosHEPEVTParticle::setPy(double py){
  m_py = py;
}


void PhotosHEPEVTParticle::setPz(double pz){
  m_pz = pz;
}

void PhotosHEPEVTParticle::setE(double e){
  m_e = e;
}

int PhotosHEPEVTParticle::getBarcode(){
  return m_barcode;
}

void PhotosHEPEVTParticle::setBarcode(int barcode){
  m_barcode = barcode;
}

void PhotosHEPEVTParticle::setEvent(PhotosHEPEVTEvent *event){
  m_event = event;
}

int PhotosHEPEVTParticle::getFirstMotherIndex(){
  return m_mother_start;
}

int PhotosHEPEVTParticle::getLastMotherIndex(){
  return m_mother_end;
}

int PhotosHEPEVTParticle::getFirstDaughterIndex(){
  return m_daughter_start;
}

int PhotosHEPEVTParticle::getLastDaughterIndex(){
  return m_daughter_end;
}

