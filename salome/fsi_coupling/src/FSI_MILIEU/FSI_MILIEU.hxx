#ifndef _FSI_MILIEU_HXX_
#define _FSI_MILIEU_HXX_

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2018 EDF S.A.

  This program is free software; you can redistribute it and/or modify it under
  the terms of the GNU General Public License as published by the Free Software
  Foundation; either version 2 of the License, or (at your option) any later
  version.

  This program is distributed in the hope that it will be useful, but WITHOUT
  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
  FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
  details.

  You should have received a copy of the GNU General Public License along with
  this program; if not, write to the Free Software Foundation, Inc., 51 Franklin
  Street, Fifth Floor, Boston, MA 02110-1301, USA.
*/

/*----------------------------------------------------------------------------*/

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

