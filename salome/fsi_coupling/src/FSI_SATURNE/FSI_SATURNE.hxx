
#ifndef _FSI_SATURNE_HXX_
#define _FSI_SATURNE_HXX_

#include <SALOME_Component.hh>
#include "Superv_Component_i.hxx"
#include "FSI.hh"

//COMPODEFS

//ENDDEF

class FSI_SATURNE_i: public virtual POA_FSI_ORB::FSI_SATURNE,
                       public virtual Superv_Component_i
{
  public:
    FSI_SATURNE_i(CORBA::ORB_ptr orb, PortableServer::POA_ptr poa,
              PortableServer::ObjectId * contId,
              const char *instanceName, const char *interfaceName);
    FSI_SATURNE_i(CORBA::ORB_ptr orb, PortableServer::POA_ptr poa,
              Engines::Container_ptr container,
              const char *instanceName, const char *interfaceName);
    virtual ~FSI_SATURNE_i();
    void destroy();
    CORBA::Boolean init_service(const char * service_name);
    void run(const char* app_name,CORBA::Long verbosity,CORBA::Long& retval);
};

extern "C"
{
    PortableServer::ObjectId * FSI_SATURNEEngine_factory( CORBA::ORB_ptr orb,
                                                      PortableServer::POA_ptr poa,
                                                      PortableServer::ObjectId * contId,
                                                      const char *instanceName,
                                                      const char *interfaceName);
    void yacsinit();
}
#endif

