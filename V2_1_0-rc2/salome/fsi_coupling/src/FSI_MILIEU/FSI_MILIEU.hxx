
#ifndef _FSI_MILIEU_HXX_
#define _FSI_MILIEU_HXX_

#include <SALOME_Component.hh>
#include "Superv_Component_i.hxx"
#include "FSI.hh"

//COMPODEFS

//ENDDEF

class FSI_MILIEU_i: public virtual POA_FSI_ORB::FSI_MILIEU,
                       public virtual Superv_Component_i
{
  public:
    FSI_MILIEU_i(CORBA::ORB_ptr orb, PortableServer::POA_ptr poa,
              PortableServer::ObjectId * contId,
              const char *instanceName, const char *interfaceName);
    FSI_MILIEU_i(CORBA::ORB_ptr orb, PortableServer::POA_ptr poa,
              Engines::Container_ptr container,
              const char *instanceName, const char *interfaceName);
    virtual ~FSI_MILIEU_i();
    void destroy();
    CORBA::Boolean init_service(const char * service_name);
    void inter_run(CORBA::Long NBPDTM,CORBA::Long NBSSIT,CORBA::Long ISYNCP,CORBA::Long NTCHR,CORBA::Double DTREF,CORBA::Double TTINIT,CORBA::Double EPSILO);
};

extern "C"
{
    PortableServer::ObjectId * FSI_MILIEUEngine_factory( CORBA::ORB_ptr orb,
                                                      PortableServer::POA_ptr poa,
                                                      PortableServer::ObjectId * contId,
                                                      const char *instanceName,
                                                      const char *interfaceName);
    void yacsinit();
}
#endif

