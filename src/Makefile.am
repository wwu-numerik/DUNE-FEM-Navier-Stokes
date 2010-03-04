# $Id: duneproject 5104 2008-03-06 16:24:04Z robertk $


# possible options
#LDADD = $(UG_LDFLAGS) $(AMIRAMESH_LDFLAGS) $(UG_LIBS) $(AMIRAMESH_LIBS)
#AM_CPPFLAGS = $(UG_CPPFLAGS) $(AMIRAMESH_CPPFLAGS)

noinst_PROGRAMS = dune_navier_stokes

dune_navier_stokes_SOURCES = dune_navier_stokes.cc

dune_navier_stokes_CXXFLAGS = $(MPI_CPPFLAGS) $(UG_CPPFLAGS) $(AMIRAMESH_CPPFLAGS) $(ALBERTA_CPPFLAGS)
dune_navier_stokes_LDADD = $(MPI_LDFLAGS) $(ALBERTA_LDFLAGS) $(ALBERTA_LIBS) $(AMIRAMESH_LDFLAGS) $(AMIRAMESH_LIBS) $(UG_LDFLAGS) $(UG_LIBS) $(MPI_LDFLAGS) $(DUNE_LDFLAGS) $(DUNE_LIBS)

# don't follow the full GNU-standard
# we need automake 1.5
AUTOMAKE_OPTIONS = foreign 1.5
# pass most important options when "make distcheck" is used
DISTCHECK_CONFIGURE_FLAGS = --with-dune=$(DUNEROOT) CXX="$(CXX)" CC="$(CC)"
