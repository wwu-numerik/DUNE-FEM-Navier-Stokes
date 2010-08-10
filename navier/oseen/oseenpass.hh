/** \file stokespass.hh
    \brief  stokespass.hh
 **/

#ifndef OSEENPASS_HH
#define OSEENPASS_HH

#include <dune/fem/pass/pass.hh>
#include <dune/fem/operator/matrix/spmatrix.hh>
#include <dune/fem/space/dgspace.hh>
#include <dune/fem/quadrature/caching/twistutility.hh>
#include <dune/fem/misc/l2norm.hh>

#include <dune/navier/oseen/inverse_operator.hh>

#include <dune/common/stdstreams.hh>
#include <dune/fem/operator/matrix/spmatrix.hh>

#ifndef NLOG // if we want logging, should be removed in the end
    #include <dune/stuff/printing.hh>
    #include <dune/stuff/misc.hh>
    #include <dune/stuff/logging.hh>
#endif

#include <dune/stuff/grid.hh>
#include <dune/stuff/functions.hh>

#include <dune/stuff/profiler.hh>

namespace Dune
{
	namespace NavierStokes {
		namespace NonlinearStep {

			template <class RowSpaceImp, class ColSpaceImp = RowSpaceImp>
			struct MatrixTraits : public Dune::SparseRowMatrixTraits<RowSpaceImp,ColSpaceImp> {
				struct StencilType {
					template <class T>
					static int nonZerosEstimate( T& rangeSpace_ ) {
						return 100;
					}
				};
			};

/**
 *  \brief  OseenPass
 *
 *  \todo   doc
 **/
			template <  class DiscreteModelImp,
						class PreviousPassImp,
						int PassID = 0 >
			class OseenPass
				: public Pass < DiscreteModelImp, PreviousPassImp, PassID >
			{
				public:
					//! own type
					typedef OseenPass< DiscreteModelImp, PreviousPassImp, PassID >
						ThisType;

					//! base type
					typedef Pass < DiscreteModelImp, PreviousPassImp, PassID >
						BaseType;

					//! previous pass type
					typedef PreviousPassImp
						PreviousPassType;

					//! discrete model type
					typedef DiscreteModelImp
						DiscreteModelType;

					//! volume quadrature type
					typedef typename DiscreteModelType::VolumeQuadratureType
						VolumeQuadratureType;

					//! face quadrature type
					typedef typename DiscreteModelType::FaceQuadratureType
						FaceQuadratureType;

					//! type of discrete function space wrapper
					typedef typename DiscreteModelType::DiscreteStokesFunctionSpaceWrapperType
						DiscreteStokesFunctionSpaceWrapperType;

					//! discrete function wrapper type
					typedef typename DiscreteModelType::DiscreteStokesFunctionWrapperType
						DiscreteStokesFunctionWrapperType;

					//! discrete function type for the velocity
					typedef typename DiscreteStokesFunctionWrapperType::DiscreteVelocityFunctionType
						DiscreteVelocityFunctionType;

					//! discrete function space type for the velocity
					typedef typename DiscreteVelocityFunctionType::DiscreteFunctionSpaceType
						DiscreteVelocityFunctionSpaceType;

					//! discrete function type for sigma
					typedef typename DiscreteModelType::DiscreteSigmaFunctionType
						DiscreteSigmaFunctionType;

					//! discrete function space type for sigma
					typedef typename DiscreteSigmaFunctionType::DiscreteFunctionSpaceType
						DiscreteSigmaFunctionSpaceType;

					//! discrete fucntion type for the pressure
					typedef typename DiscreteStokesFunctionWrapperType::DiscretePressureFunctionType
						DiscretePressureFunctionType;

					//! discrete function space type for the pressure
					typedef typename DiscretePressureFunctionType::DiscreteFunctionSpaceType
						DiscretePressureFunctionSpaceType;

					typedef typename DiscreteModelType::Traits::ExtraDataFunctionType
						ExtraDataFunctionType;

					//! Coordinate type on the element
					typedef typename DiscreteVelocityFunctionSpaceType::DomainType
						ElementCoordinateType;

					//! Coordinate type on an intersection
					typedef typename FaceQuadratureType::LocalCoordinateType
						IntersectionCoordinateType;

					//! Vector type of the velocity's discrete function space's range
					typedef typename DiscreteVelocityFunctionSpaceType::RangeType
						VelocityRangeType;

					typedef typename DiscreteVelocityFunctionSpaceType::BaseFunctionSetType::JacobianRangeType
						VelocityJacobianRangeType;

					//! vector type of sigmas' discrete functions space's range
					typedef typename DiscreteSigmaFunctionSpaceType::RangeType
						SigmaRangeType;

					typedef typename DiscreteSigmaFunctionSpaceType::BaseFunctionSetType::JacobianRangeType
						SigmaJacobianRangeType;

					//! Vector type of the pressure's discrete function space's range
					typedef typename DiscretePressureFunctionSpaceType::RangeType
						PressureRangeType;

					typedef typename DiscretePressureFunctionSpaceType::BaseFunctionSetType::JacobianRangeType
						PressureJacobianRangeType;

					//! Type of GridPart
					typedef typename DiscreteVelocityFunctionSpaceType::GridPartType
						GridPartType;

					//! Intersection iterator of the gridpart
					typedef typename GridPartType::IntersectionIteratorType
						IntersectionIteratorType;

					//! local coordinate type on an intersection
					typedef typename FaceQuadratureType::LocalCoordinateType
						LocalIntersectionCoordinateType;

					//! entity iterator of the gridpart
					typedef typename GridPartType::template Codim< 0 >::IteratorType
						EntityIteratorType;

					//! type of the grid
					typedef typename GridPartType::GridType
						GridType;

					//! type of codim 0 entity
					typedef typename GridType::template Codim< 0 >::Entity
						EntityType;

			#ifdef USE_NESTED_CG_SOLVER
					//! alternative, not fully functional, solver implemementation
					typedef NestedCgSaddlepointInverseOperator< ThisType >
						AltInvOpType;
			#endif
					//! type of the used solver
					typedef ReducedInverseOperator< ThisType >
						InvOpType;


					//! polynomial order for the discrete sigma function space
					static const int sigmaSpaceOrder
						= DiscreteModelType::sigmaSpaceOrder;
					//! polynomial order for the discrete velocity function space
					static const int velocitySpaceOrder
						= DiscreteModelType::velocitySpaceOrder;
					//! polynomial order for the discrete pressure function space
					static const int pressureSpaceOrder
						= DiscreteModelType::pressureSpaceOrder;


					/**
					 *  \name typedefs needed for interface compliance
					 *  \{
					 **/
					typedef typename BaseType::DestinationType
						DestinationType;

					typedef typename BaseType::DomainType
						DomainType;

					typedef typename BaseType::RangeType
						RangeType;

					typedef typename BaseType::TotalArgumentType
						TotalArgumentType;
					/**
					 *  \}
					 **/

					/**
					 *  \brief  constructor
					 *  \todo   doc
					 **/
					OseenPass( PreviousPassType& prevPass,
								DiscreteModelType& discreteModel,
								GridPartType& gridPart,
								DiscreteStokesFunctionSpaceWrapperType& spaceWrapper,
								const ExtraDataFunctionType& beta )
						: BaseType( prevPass ),
						discreteModel_( discreteModel ),
						gridPart_( gridPart ),
						spaceWrapper_( spaceWrapper ),
						velocitySpace_( spaceWrapper.discreteVelocitySpace() ),
						pressureSpace_( spaceWrapper.discretePressureSpace() ),
						sigmaSpace_( gridPart ),
						beta_( beta )
					{}

					/**
					 *  \brief  empty constructor
					 **/
					OseenPass()
					{}

					//! used in Postprocessing to get refs to gridparts, spaces
					const DiscreteStokesFunctionSpaceWrapperType& GetFunctionSpaceWrapper() const
					{
						return spaceWrapper_;
					}

					/**
					 *  \todo doc
					 *  \attention  think about quadrature orders
					 **/
					virtual void apply( const DomainType &arg, RangeType &dest) const
					{

						// profiler information
						profiler().StartTiming("Pass -- ASSEMBLE");

						// entity and geometry types
						typedef typename EntityType::Geometry
							EntityGeometryType;
						typedef typename Dune::FieldMatrix< typename EntityGeometryType::ctype,
															EntityGeometryType::coorddimension,
															EntityGeometryType::mydimension >
							JacobianInverseTransposedType;

						// viscosity
						const double mu = discreteModel_.viscosity();

						// generalized stokes alpha
						const double alpha = discreteModel_.alpha();

						// matrices
						// M\in R^{M\times M}
						typedef SparseRowMatrixObject<  DiscreteSigmaFunctionSpaceType,
														DiscreteSigmaFunctionSpaceType,
														MatrixTraits<DiscreteSigmaFunctionSpaceType,DiscreteSigmaFunctionSpaceType> >
							MInversMatrixType;
						MInversMatrixType MInversMatrix( sigmaSpace_, sigmaSpace_ );
						MInversMatrix.reserve();
						assert( MInversMatrix.matrix().rows() == MInversMatrix.matrix().cols() );
						// W\in R^{M\times L}
						typedef SparseRowMatrixObject<  DiscreteSigmaFunctionSpaceType,
														DiscreteVelocityFunctionSpaceType,
														MatrixTraits<DiscreteSigmaFunctionSpaceType, DiscreteVelocityFunctionSpaceType> >
							WmatrixType;
						WmatrixType Wmatrix( sigmaSpace_, velocitySpace_ );
						Wmatrix.reserve();
						// X\in R^{L\times M}
						typedef SparseRowMatrixObject<  DiscreteVelocityFunctionSpaceType,
														DiscreteSigmaFunctionSpaceType,
														MatrixTraits<DiscreteVelocityFunctionSpaceType, DiscreteSigmaFunctionSpaceType> >
							XmatrixType;
						XmatrixType Xmatrix( velocitySpace_, sigmaSpace_ );
						Xmatrix.reserve();
						// Y\in R^{L\times L}
						typedef SparseRowMatrixObject<  DiscreteVelocityFunctionSpaceType,
														DiscreteVelocityFunctionSpaceType,
														MatrixTraits<DiscreteVelocityFunctionSpaceType,DiscreteVelocityFunctionSpaceType> >
							YmatrixType;
						YmatrixType Ymatrix( velocitySpace_, velocitySpace_ );
						Ymatrix.reserve();
						// Z\in R^{L\times K}
						typedef SparseRowMatrixObject<  DiscreteVelocityFunctionSpaceType,
														DiscretePressureFunctionSpaceType,
														MatrixTraits<DiscreteVelocityFunctionSpaceType,DiscretePressureFunctionSpaceType> >
							ZmatrixType;
						ZmatrixType Zmatrix( velocitySpace_, pressureSpace_ );
						Zmatrix.reserve();
						// E\in R^{K\times L}
						typedef SparseRowMatrixObject<  DiscretePressureFunctionSpaceType,
														DiscreteVelocityFunctionSpaceType,
														MatrixTraits<DiscretePressureFunctionSpaceType,DiscreteVelocityFunctionSpaceType> >
							EmatrixType;
						EmatrixType Ematrix( pressureSpace_, velocitySpace_ );
						Ematrix.reserve();
						// R\in R^{K\times K}
						typedef SparseRowMatrixObject<  DiscretePressureFunctionSpaceType,
														DiscretePressureFunctionSpaceType,
														MatrixTraits<DiscretePressureFunctionSpaceType,DiscretePressureFunctionSpaceType> >
							RmatrixType;
						RmatrixType Rmatrix( pressureSpace_, pressureSpace_ );
						Rmatrix.reserve();

						// local matrices
						// M\in R^{M\times M}
						typedef typename MInversMatrixType::LocalMatrixType
							LocalMInversMatrixType;
						// W\in R^{M\times L}
						typedef typename WmatrixType::LocalMatrixType
							LocalWmatrixType;
						// X\in R^{L\times M}
						typedef typename XmatrixType::LocalMatrixType
							LocalXmatrixType;
						// Y\in R^{L\times L}
						typedef typename YmatrixType::LocalMatrixType
							LocalYmatrixType;
						// Z\in R^{L\times K}
						typedef typename ZmatrixType::LocalMatrixType
							LocalZmatrixType;
						// E\in R^{K\times L}
						typedef typename EmatrixType::LocalMatrixType
							LocalEmatrixType;
						// R\in R^{K\times K}
						typedef typename RmatrixType::LocalMatrixType
							LocalRmatrixType;

						// right hand sides
						// H_{1}\in R^{M}
						DiscreteSigmaFunctionType H1rhs( "H1", sigmaSpace_ );
						H1rhs.clear();
						// H_{2}\in R^{L}
						DiscreteVelocityFunctionType H2rhs( "H2", velocitySpace_ );
						H2rhs.clear();
						// H_{3}\in R^{K}
						DiscretePressureFunctionType H3rhs( "H3", pressureSpace_ );
						H3rhs.clear();

						// local right hand sides
						// H_{1}\in R^{M}
						typedef typename DiscreteSigmaFunctionType::LocalFunctionType
							LocalH1rhsType;
						// H_{2}\in R^{L}
						typedef typename DiscreteVelocityFunctionType::LocalFunctionType
							LocalH2rhsType;
						// H_{3}\in R^{K}
						typedef typename DiscretePressureFunctionType::LocalFunctionType
							LocalH3rhsType;

						// base functions
						// of type sigma
						typedef typename DiscreteSigmaFunctionSpaceType::BaseFunctionSetType
							SigmaBaseFunctionSetType;
						// of type u
						typedef typename DiscreteVelocityFunctionSpaceType::BaseFunctionSetType
							VelocityBaseFunctionSetType;
						// of type p
						typedef typename DiscretePressureFunctionSpaceType::BaseFunctionSetType
							PressureBaseFunctionSetType;

						// eps
						const double eps = Parameters().getParam( "eps", 1.0e-14 );

			#ifndef NLOG
						// logging stuff
						Logging::LogStream& infoStream = Logger().Info();
			//            Logging::LogStream& debugStream = Logger().Info();
						Logging::LogStream& debugStream = Logger().Dbg(); // sometimes Dbg() doesn't work
						bool entityOutput = false;
						bool intersectionOutput = false;
						const int outputEntity = 0;
						const int outputIntersection = -1;
						int entityNR = 0;
						int intersectionNR = 0;
						int numberOfEntities = 0;
						int numberOfIntersections = 0;
						int numberOfBoundaryIntersections = 0;
						int numberOfInnerIntersections = 0;
						const bool Mprint = Parameters().getParam( "Mprint", false );
						const bool Wprint = Parameters().getParam( "Wprint", false );
						const bool Xprint = Parameters().getParam( "Xprint", false );
						const bool Yprint = Parameters().getParam( "Yprint", false );
						const bool Zprint = Parameters().getParam( "Zprint", false );
						const bool Eprint = Parameters().getParam( "Eprint", false );
						const bool Rprint = Parameters().getParam( "Rprint", false );
						const bool H1print = Parameters().getParam( "H1print", false );
						const bool H2print = Parameters().getParam( "H2print", false );
						const bool H3print = Parameters().getParam( "H3print", false );
						const bool allOutput = Parameters().getParam( "allOutput", false );
						const bool Mdebug = Parameters().getParam( "Mdebug", false );
						const bool Wdebug = Parameters().getParam( "Wdebug", false );
						const bool Xdebug = Parameters().getParam( "Xdebug", false );
						const bool Ydebug = Parameters().getParam( "Ydebug", false );
						const bool Zdebug = Parameters().getParam( "Zdebug", false );
						const bool Edebug = Parameters().getParam( "Edebug", false );
						const bool Rdebug = Parameters().getParam( "Rdebug", false );
						const bool H1debug = Parameters().getParam( "H1debug", false );
						const bool H2debug = Parameters().getParam( "H2debug", false );
						const bool H3debug = Parameters().getParam( "H3debug", false );
						int fivePercentOfEntities = 0;
						int fivePercents = 0;
						infoStream << "this is OseenPass::apply()" << std::endl;

						// do an empty grid walk to get informations
						double maxGridWidth( 0.0 );
						EntityIteratorType entityItEndLog = velocitySpace_.end();
						for (   EntityIteratorType entityItLog = velocitySpace_.begin();
								entityItLog != entityItEndLog;
								++entityItLog ) {
							const EntityType& entity = *entityItLog;
							// count entities
							++numberOfEntities;
							// walk the intersections
							IntersectionIteratorType intItEnd = gridPart_.iend( entity );
							for (   IntersectionIteratorType intIt = gridPart_.ibegin( entity );
									intIt != intItEnd;
									++intIt ) {
								// count intersections
								++numberOfIntersections;
								maxGridWidth = std::max( Stuff::getLenghtOfIntersection( *intIt ), maxGridWidth );
								// if we are inside the grid
								if ( intIt->neighbor() && !intIt->boundary() ) {
									// count inner intersections
									++numberOfInnerIntersections;
								}
								// if we are on the boundary of the grid
								if ( !intIt->neighbor() && intIt->boundary() ) {
									// count boundary intersections
									++numberOfBoundaryIntersections;
								}
							}
						}
						if ( numberOfEntities > 19 ) {
							infoStream << "found " << numberOfEntities << " entities," << std::endl;
							infoStream << "found " << numberOfIntersections << " intersections," << std::endl;
							infoStream << "      " << numberOfInnerIntersections << " intersections inside and" << std::endl;
							infoStream << "      " << numberOfBoundaryIntersections << " intersections on the boundary." << std::endl;
							infoStream << "      maxGridWidth is " << maxGridWidth << std::endl;
							infoStream << "- starting gridwalk" << std::endl;
							fivePercentOfEntities = int( std::floor(double(numberOfEntities) / double(20)));
							infoStream << "  [ assembling         ]" << std::endl;
							infoStream << "  [";
						} else {
							infoStream << "found " << numberOfEntities << " entities," << std::endl;
							infoStream << "found " << numberOfIntersections << " intersections," << std::endl;
							infoStream << "      " << numberOfInnerIntersections << " intersections inside and" << std::endl;
							infoStream << "      " << numberOfBoundaryIntersections << " intersections on the boundary." << std::endl;
							infoStream << "      maxGridWidth is " << maxGridWidth << std::endl;
							infoStream << "- starting gridwalk" << std::endl;
						}
						infoStream.Suspend();
			#endif

						// walk the grid
						EntityIteratorType entityItEnd = velocitySpace_.end();
						for (   EntityIteratorType entityIt = velocitySpace_.begin();
								entityIt != entityItEnd;
								++entityIt ) {

							// get entity and geometry
							const EntityType& entity = *entityIt;
							const EntityGeometryType& geometry = entity.geometry();

							// get local matrices for the volume integral
							LocalMInversMatrixType localMInversMatrixElement = MInversMatrix.localMatrix( entity, entity );
							LocalWmatrixType localWmatrixElement = Wmatrix.localMatrix( entity, entity );
							LocalXmatrixType localXmatrixElement = Xmatrix.localMatrix( entity, entity );
							LocalYmatrixType localYmatrixElement = Ymatrix.localMatrix( entity, entity );
							LocalZmatrixType localZmatrixElement = Zmatrix.localMatrix( entity, entity );
							LocalEmatrixType localEmatrixElement = Ematrix.localMatrix( entity, entity );
							LocalRmatrixType localRmatrixElement = Rmatrix.localMatrix( entity, entity );

							// get local right hand sides
							LocalH1rhsType localH1rhs = H1rhs.localFunction( entity );
							LocalH2rhsType localH2rhs = H2rhs.localFunction( entity );
							LocalH3rhsType localH3rhs = H3rhs.localFunction( entity );

							// get basefunctionsets
							const SigmaBaseFunctionSetType sigmaBaseFunctionSetElement = sigmaSpace_.baseFunctionSet( entity );
							const VelocityBaseFunctionSetType velocityBaseFunctionSetElement = velocitySpace_.baseFunctionSet( entity );
							const PressureBaseFunctionSetType pressureBaseFunctionSetElement = pressureSpace_.baseFunctionSet( entity );
							const int numSigmaBaseFunctionsElement = sigmaBaseFunctionSetElement.numBaseFunctions();
							const int numVelocityBaseFunctionsElement = velocityBaseFunctionSetElement.numBaseFunctions();
							const int numPressureBaseFunctionsElement = pressureBaseFunctionSetElement.numBaseFunctions();

							// get quadrature
							const VolumeQuadratureType volumeQuadratureElement( entity,
																				( 4 * pressureSpaceOrder ) + 1 );

			#ifndef NLOG
							if ( numberOfEntities > 19 ) {
								if ( ( entityNR % fivePercentOfEntities ) == 0 ) {
									if ( fivePercents < 20 ) {
										infoStream.Resume();
										infoStream << "=";
										infoStream.Flush();
										infoStream.Suspend();
										++fivePercents;
									}
								}
							}
							debugStream.Suspend(); // disable logging
			//                if ( outputEntity == entityNR ) entityOutput = true;
							if ( allOutput ) entityOutput = true;
							if ( entityOutput ) debugStream.Resume(); // enable logging
							debugStream << "  - numSigmaBaseFunctionsElement: " << numSigmaBaseFunctionsElement << std::endl;
							debugStream << "  - numVelocityBaseFunctionsElement: " << numVelocityBaseFunctionsElement << std::endl;
							debugStream << "  - numPressureBaseFunctionsElement: " << numPressureBaseFunctionsElement << std::endl;
							debugStream << "  - == start calculations on entity " << entityNR << std::endl;
							debugStream << "    - entity " << entityNR << " has " << geometry.corners() << " corners:";
							for ( int i = 0; i < geometry.corners(); ++i ) {
								const VelocityRangeType corner = geometry.corner( i );
								Stuff::printFieldVector( corner, "corner_"+Stuff::toString(i), debugStream, "      " );
							}
							debugStream << std::endl;
							bool Moutput = false;
							bool Woutput = false;
							bool Xoutput = false;
							bool Youtput = false;
							bool Zoutput = false;
							bool Eoutput = false;
							bool Routput = false;
							bool H1output = false;
							bool H2output = false;
							bool H3output = false;
							// we want logging at the following base functions
							const int logBaseI = Parameters().getParam( "logBaseI", 0 );
							const int logBaseJ = Parameters().getParam( "logBaseJ", 0 );
							debugStream.Suspend(); // disable logging
			#endif

							// compute volume integrals

							//                                                     // we will call this one
							// (M^{-1})_{i,j} = (\int_{T}\tau_{j}:\tau_{i}dx)^{-1} // Minvs' volume integral
			#ifndef NLOG
							if ( Mdebug ) {
								debugStream.Resume(); // enable logging
								debugStream << "    = M volume =================" << std::endl;
								debugStream.Suspend();
							}
			#endif
							for ( int i = 0; i < numSigmaBaseFunctionsElement; ++i ) {
								for ( int j = 0; j < numSigmaBaseFunctionsElement; ++j ) {
									double M_i_j = 0.0;
			#ifndef NLOG
			//                        if ( ( i == logBaseI ) && ( j == logBaseJ ) ) Moutput = true;
			//                        if ( allOutput ) Moutput = true;
									if ( Mdebug ) Moutput = true;
									if ( entityOutput && Moutput ) debugStream.Resume(); // enable logging
									debugStream << "    = M ========================" << std::endl;
									debugStream << "    basefunctions " << i << " " << j << std::endl;
									debugStream << "    volumeQuadrature.nop() " << volumeQuadratureElement.nop() << std::endl;
			#endif
									// sum over all quadrature points
									for ( size_t quad = 0; quad < volumeQuadratureElement.nop(); ++quad ) {
										// get x
										const ElementCoordinateType x = volumeQuadratureElement.point( quad );
										// get the integration factor
										const double elementVolume = geometry.integrationElement( x );
										// get the quadrature weight
										const double integrationWeight = volumeQuadratureElement.weight( quad );
										// compute \tau_{i}:\tau_{j}
										SigmaRangeType tau_i( 0.0 );
										SigmaRangeType tau_j( 0.0 );
										sigmaBaseFunctionSetElement.evaluate( i, x, tau_i );
										sigmaBaseFunctionSetElement.evaluate( j, x, tau_j );
										const double tau_i_times_tau_j = Stuff::colonProduct( tau_i, tau_j );
										// compute M_i_j
										M_i_j += elementVolume
											* integrationWeight
											* tau_i_times_tau_j;
			#ifndef NLOG
										debugStream << "    - quadPoint " << quad;
										Stuff::printFieldVector( x, "x", debugStream, "      " );
										debugStream << "\n        - elementVolume: " << elementVolume << std::endl;
										debugStream << "        - integrationWeight: " << integrationWeight;
										Stuff::printFieldMatrix( tau_i, "tau_i", debugStream, "      " );
										Stuff::printFieldMatrix( tau_j, "tau_j", debugStream, "      " );
										debugStream << "\n        - tau_i_times_tau_j: " << tau_i_times_tau_j << std::endl;
										debugStream << "        - M_" << i << "_" << j << "+=: " << M_i_j << std::endl;
			#endif
									} // done sum over quadrature points
									// if small, should be zero
									if ( fabs( M_i_j ) < eps ) {
										M_i_j = 0.0;
									} // else invert
									else {
										M_i_j = 1.0 / M_i_j;
										// add to matrix
										localMInversMatrixElement.add( i, j, M_i_j );
									}
			#ifndef NLOG
									Moutput = false;
									debugStream.Suspend(); // disable logging
			#endif
								}
							} // done computing Minvs' volume integral

							//                                                        // we will call this one
							// (W)_{i,j} += \int_{T}v_{j}\cdot(\nabla\cdot\tau_{i})dx // W's volume integral
							//                                                        // see also "W's entitity surface integral", "W's neighbour surface integral" and "W's boundary integral" below
			#ifndef NLOG
							if ( Wdebug ) {
								debugStream.Resume(); // enable logging
								debugStream << "    = W volume =================" << std::endl;
								debugStream.Suspend();
							}
			#endif
							for ( int i = 0; i < numSigmaBaseFunctionsElement; ++i ) {
								for ( int j = 0; j < numVelocityBaseFunctionsElement; ++j ) {
									double W_i_j = 0.0;
			#ifndef NLOG
			//                        if ( ( i == logBaseI ) && ( j == logBaseJ ) ) Woutput = true;
			//                        if ( ( ( i == logBaseI ) && ( j == logBaseJ ) ) && Wdebug ) Woutput = true;
			//                        if ( allOutput ) Woutput = true;
									if ( Wdebug ) Woutput = true;
									if ( entityOutput && Woutput ) debugStream.Resume(); // enable logging
									debugStream << "    = W ========================" << std::endl;
									debugStream << "    basefunctions " << i << " " << j << std::endl;
									debugStream << "    volumeQuadrature.nop() " << volumeQuadratureElement.nop() << std::endl;
			#endif
									// sum over all quadrature points
									for ( size_t quad = 0; quad < volumeQuadratureElement.nop(); ++quad ) {
										// get x
										const ElementCoordinateType x = volumeQuadratureElement.point( quad );
										// get the integration factor
										const double elementVolume = geometry.integrationElement( x );
										// get the quadrature weight
										const double integrationWeight = volumeQuadratureElement.weight( quad );
										// compute v_j^t \cdot ( \nabla \cdot \tau_i^t )
										VelocityRangeType v_j( 0.0 );
										velocityBaseFunctionSetElement.evaluate( j, x, v_j );
										const double divergence_of_tau_i_times_v_j = sigmaBaseFunctionSetElement.evaluateGradientSingle( i, entity, x, prepareVelocityRangeTypeForSigmaDivergence( v_j ) );
										W_i_j += elementVolume
											* integrationWeight
											* divergence_of_tau_i_times_v_j;
			#ifndef NLOG
										debugStream << "    - quadPoint " << quad;
										Stuff::printFieldVector( x, "x", debugStream, "      " );
										debugStream << "\n        - elementVolume: " << elementVolume << std::endl;
										debugStream << "        - integrationWeight: " << integrationWeight;
										Stuff::printFieldVector( v_j, "v_j", debugStream, "      " );
										debugStream << "\n        - divergence_of_tau_i_times_v_j: " << divergence_of_tau_i_times_v_j << std::endl;
										debugStream << "        - W_" << i << "_" << j << "+=: " << W_i_j << std::endl;
			#endif
									} // done sum over quadrature points
									// if small, should be zero
									if ( fabs( W_i_j ) < eps ) {
										W_i_j = 0.0;
									}
									else
										// add to matrix
										localWmatrixElement.add( i, j, W_i_j );
			#ifndef NLOG
									Woutput = false;
									debugStream.Suspend(); // disable logging
			#endif
								}
							} // done computing W's volume integral

							//                                                  // we will call this one
							// (X)_{i,j} += \mu\int_{T}\tau_{j}:\nabla v_{i} dx // X's volume integral
							//                                                  // see also "X's entitity surface integral", "X's neighbour surface integral" and "X's boundary integral" below
			#ifndef NLOG
							if ( Xdebug ) {
								debugStream.Resume(); // enable logging
								debugStream << "    = X volume =================" << std::endl;
								debugStream.Suspend();
							}
			#endif
							for ( int i = 0; i < numVelocityBaseFunctionsElement; ++i ) {
								for ( int j = 0; j < numSigmaBaseFunctionsElement; ++j ) {
									double X_i_j = 0.0;
			#ifndef NLOG
			//                        if ( ( i == logBaseI ) && ( j == logBaseJ ) ) Xoutput = true;
			//                        if ( allOutput ) Xoutput = true;
									if ( Xdebug ) Xoutput = true;
									if ( entityOutput && Xoutput ) debugStream.Resume(); // enable logging
									debugStream << "    = X ========================" << std::endl;
									debugStream << "    basefunctions " << i << " " << j << std::endl;
									debugStream << "    volumeQuadrature.nop() " << volumeQuadratureElement.nop() << std::endl;
			#endif
									// sum over all quadrature points
									for ( size_t quad = 0; quad < volumeQuadratureElement.nop(); ++quad ) {
										// get x
										const ElementCoordinateType x = volumeQuadratureElement.point( quad );
										// get the integration factor
										const double elementVolume = geometry.integrationElement( x );
										// get the quadrature weight
										const double integrationWeight = volumeQuadratureElement.weight( quad );
										// compute \tau_{j}:\nabla v_{i}
										SigmaRangeType tau_j( 0.0 );
										sigmaBaseFunctionSetElement.evaluate( j, x, tau_j );
										const double gradient_of_v_i_times_tau_j = velocityBaseFunctionSetElement.evaluateGradientSingle( i, entity, x, tau_j );
										X_i_j += elementVolume
											* integrationWeight
											* mu
											* gradient_of_v_i_times_tau_j;
			#ifndef NLOG
										debugStream << "    - quadPoint " << quad;
										Stuff::printFieldVector( x, "x", debugStream, "      " );
										debugStream << "\n        - elementVolume: " << elementVolume << std::endl;
										debugStream << "        - integrationWeight: " << integrationWeight;
										Stuff::printFieldMatrix( tau_j, "tau_j", debugStream, "      " );
										debugStream << "\n        - gradient_of_v_i_times_tau_j: " << gradient_of_v_i_times_tau_j << std::endl;
										debugStream << "        - X_" << i << "_" << j << "+=: " << X_i_j << std::endl;
			#endif
									} // done sum over quadrature points
									// if small, should be zero
									if ( fabs( X_i_j ) < eps ) {
										X_i_j = 0.0;
									}
									else
										// add to matrix
										localXmatrixElement.add( i, j, X_i_j );
			#ifndef NLOG
									Xoutput = false;
									debugStream.Suspend(); // disable logging
			#endif
								}
							} // done computing X's volume integral

							//
							// (Y)
							//
			#ifndef NLOG
							if ( Ydebug ) {
								debugStream.Resume(); // enable logging
								debugStream << "    = Y volume =================" << std::endl;
								debugStream.Suspend();
							}
			#endif
			//                if ( discreteModel_.isGeneralized() )
							{
							for ( int i = 0; i < numVelocityBaseFunctionsElement; ++i ) {
								for ( int j = 0; j < numVelocityBaseFunctionsElement; ++j ) {
									double Y_i_j = 0.0;
			#ifndef NLOG
			//                        if ( ( i == logBaseI ) && ( j == logBaseJ ) ) Xoutput = true;
			//                        if ( allOutput ) Youtput = true;
									if ( Ydebug ) Youtput = true;
									if ( entityOutput && Youtput ) debugStream.Resume(); // enable logging
									debugStream << "    = Y ========================" << std::endl;
									debugStream << "    basefunctions " << i << " " << j << std::endl;
									debugStream << "    volumeQuadrature.nop() " << volumeQuadratureElement.nop() << std::endl;
			#endif
									// sum over all quadrature points
									for ( size_t quad = 0; quad < volumeQuadratureElement.nop(); ++quad ) {
										// get x
										const ElementCoordinateType x = volumeQuadratureElement.point( quad );
										const VelocityRangeType xWorld = geometry.global( x );
										// get the integration factor
										const double elementVolume = geometry.integrationElement( x );
										// get the quadrature weight
										const double integrationWeight = volumeQuadratureElement.weight( quad );
										// compute \tau_{j}:\nabla v_{i}
										VelocityRangeType v_i( 0.0 );
										velocityBaseFunctionSetElement.evaluate( i, x, v_i );
										VelocityRangeType v_j( 0.0 );
										velocityBaseFunctionSetElement.evaluate( j, x, v_j );
										const double v_i_times_v_j = v_i * v_j;
										Y_i_j += elementVolume
											* integrationWeight
											* alpha
											* v_i_times_v_j;

										//calc u_h * \nabla * (v \tensor \beta )
										VelocityRangeType beta_eval;
										beta_.localFunction( entity ).evaluate( x, beta_eval );

//										VelocityRangeType derep = divergence_of_dyadic_product(  );
										VelocityJacobianRangeType beta_jacobian,v_i_jacobian;
//										beta_.jacobian( xWorld, beta_jacobian );
										velocityBaseFunctionSetElement.jacobian( i, x, beta_jacobian );
										const typename DiscreteVelocityFunctionType::LocalFunctionType& beta_lf =
												beta_.localFunction( entity );
										for ( size_t l = 0; l < beta_eval.dim(); ++l ) {
											for ( size_t m = 0; m < beta_eval.dim(); ++m ) {
												beta_jacobian[l][m] *= beta_lf[l];
											}
										}
//										Stuff::printFieldMatrix( beta_jacobian, "jaco", std::cerr, "JACOB ");

										velocityBaseFunctionSetElement.jacobian( i, x, v_i_jacobian );
										VelocityRangeType divergence_of_beta_v_j_tensor_beta;
										for ( size_t l = 0; l < beta_eval.dim(); ++l ) {
											double row_result = 0;
											for ( size_t m = 0; m < beta_eval.dim(); ++m ) {
												row_result += beta_jacobian[l][m] * v_i[l] + v_i_jacobian[l][m] * beta_eval[l];
											}
											divergence_of_beta_v_j_tensor_beta[l] = row_result;
										}
//										for ( size_t l = 0; l < beta_eval.dim(); ++l ) {
//											assert( !isnan(divergence_of_beta_v_j_tensor_beta[l]) );
//										}
//										Stuff::printFieldVector( divergence_of_beta_v_j_tensor_beta, "jaco", std::cerr, "DIV ");
										const double u_h_times_divergence_of_beta_v_j_tensor_beta =
												v_j * divergence_of_beta_v_j_tensor_beta;
										Y_i_j += elementVolume
											* integrationWeight
											* (-1) * u_h_times_divergence_of_beta_v_j_tensor_beta;

			#ifndef NLOG
										debugStream << "    - quadPoint " << quad;
										Stuff::printFieldVector( x, "x", debugStream, "      " );
										debugStream << "\n        - elementVolume: " << elementVolume << std::endl;
										debugStream << "        - integrationWeight: " << integrationWeight;
										Stuff::printFieldVector( v_i, "v_i", debugStream, "      " );
										Stuff::printFieldVector( v_j, "v_j", debugStream, "      " );
										debugStream << "\n        - v_i_times_v_j: " << v_i_times_v_j << std::endl;
										debugStream << "        - Y_" << i << "_" << j << "+=: " << Y_i_j << std::endl;
			#endif
									} // done sum over quadrature points
									// if small, should be zero
									if ( fabs( Y_i_j ) < eps ) {
										Y_i_j = 0.0;
									}
									else
										// add to matrix
										localYmatrixElement.add( i, j, Y_i_j );
			#ifndef NLOG
									Youtput = false;
									debugStream.Suspend(); // disable logging
			#endif
								}
							} // done computing Y's volume integral
							}

							//                                                  // we will call this one
							// (Z)_{i,j} += -\int_{T}q_{j}(\nabla\cdot v_{i})dx // Z's volume integral
							//                                                  // see also "Z's entitity surface integral", "Z's neighbour surface integral" and "Z's boundary integral" below
			#ifndef NLOG
							if ( Zdebug ) {
								debugStream.Resume(); // enable logging
								debugStream << "    = Z volume =================" << std::endl;
								debugStream.Suspend();
							}
			#endif
							for ( int i = 0; i < numVelocityBaseFunctionsElement; ++i ) {
								for ( int j = 0; j < numPressureBaseFunctionsElement; ++j ) {
									double Z_i_j = 0.0;
			#ifndef NLOG
			//                        if ( ( i == logBaseI ) && ( j == logBaseJ ) ) Zoutput = true;
			//                        if ( allOutput ) Zoutput = true;
									if ( Zdebug ) Zoutput = true;
									if ( entityOutput && Zoutput ) debugStream.Resume(); // enable logging
									debugStream << "    = Z ========================" << std::endl;
									debugStream << "    basefunctions " << i << " " << j << std::endl;
									debugStream << "    volumeQuadrature.nop() " << volumeQuadratureElement.nop() << std::endl;
			#endif
									// sum over all quadratur points
									for ( size_t quad = 0; quad < volumeQuadratureElement.nop(); ++ quad ) {
										// get x
										const ElementCoordinateType x = volumeQuadratureElement.point( quad );
										// get the integration factor
										const double elementVolume = geometry.integrationElement( x );
										// get the quadrature weight
										const double integrationWeight = volumeQuadratureElement.weight( quad );
										// compute q_{j}\cdot(\nabla\cdot v_i)
										PressureRangeType q_j( 0.0 );
										pressureBaseFunctionSetElement.evaluate( j, x, q_j );
										const double divergence_of_v_i_times_q_j = velocityBaseFunctionSetElement.evaluateGradientSingle( i, entity, x, preparePressureRangeTypeForVelocityDivergence( q_j ) );
										Z_i_j += -1.0
											* elementVolume
											* integrationWeight
											* divergence_of_v_i_times_q_j;
			#ifndef NLOG
										debugStream << "    - quadPoint " << quad;
										Stuff::printFieldVector( x, "x", debugStream, "      " );
										debugStream << "\n        - elementVolume: " << elementVolume << std::endl;
										debugStream << "        - integrationWeight: " << integrationWeight;
										Stuff::printFieldVector( q_j, "q_j", debugStream, "      " );
										debugStream << "\n        - divergence_of_v_i_times_q_j: " << divergence_of_v_i_times_q_j << std::endl;
										debugStream << "        - Z_" << i << "_" << j << "+=: " << Z_i_j << std::endl;
			#endif
									} // done sum over all quadrature points
									// if small, should be zero
									if ( fabs( Z_i_j ) < eps ) {
										Z_i_j = 0.0;
									}
									else
										// add to matrix
										localZmatrixElement.add( i, j, Z_i_j );
			#ifndef NLOG
									Zoutput = false;
									debugStream.Suspend(); // disable logging
			#endif
								}
							} // done computing Z's volume integral

							//                                    // we will call this one
							// (H2)_{j} += \int_{T}f\cdot v_{j}dx // H2's volume integral
							//                                    // see also "H2's boundary integral" further down
							if ( discreteModel_.hasForce() ) {
								for ( int j = 0; j < numVelocityBaseFunctionsElement; ++j ) {
									double H2_j = 0.0;
			#ifndef NLOG
			//                        if ( ( j == logBaseJ ) ) H2output = true;
			//                        if ( allOutput ) H2output = true;
									if ( H2debug ) H2output = true;
									if ( entityOutput && H2output ) debugStream.Resume(); // enable logging
									debugStream << "    = H2 =======================" << std::endl;
									debugStream << "    basefunction " << " " << j << std::endl;
									debugStream << "    volumeQuadrature.nop() " << volumeQuadratureElement.nop() << std::endl;
			#endif
									// sum over all quadratur points
									for ( size_t quad = 0; quad < volumeQuadratureElement.nop(); ++ quad ) {
										// get x
										const ElementCoordinateType x = volumeQuadratureElement.point( quad );
										const VelocityRangeType xWorld = geometry.global( x );
										// get the integration factor
										const double elementVolume = geometry.integrationElement( x );
										// get the quadrature weight
										const double integrationWeight = volumeQuadratureElement.weight( quad );
										// compute f\cdot v_j
										VelocityRangeType v_j( 0.0 );
										velocityBaseFunctionSetElement.evaluate( j, x, v_j );
										VelocityRangeType f( 0.0 );
#if MODEL_PROVIDES_LOCALFUNCTION
										discreteModel_.forceF().localFunction(entity).evaluate( x, f );
#else
										discreteModel_.force( 0.0, xWorld, f );
#endif
										const double f_times_v_j = f * v_j;
										H2_j += elementVolume
											* integrationWeight
											* f_times_v_j;
			#ifndef NLOG
										debugStream << "    - quadPoint " << quad;
										Stuff::printFieldVector( x, "x", debugStream, "      " );
										Stuff::printFieldVector( xWorld, "xWorld", debugStream, "      " );
										debugStream << "\n        - elementVolume: " << elementVolume << std::endl;
										debugStream << "        - integrationWeight: " << integrationWeight;
										Stuff::printFieldVector( f, "f", debugStream, "      " );
										Stuff::printFieldVector( v_j, "v_j", debugStream, "      " );
										debugStream << "\n        - f_times_v_j: " << f_times_v_j << std::endl;
										debugStream << "        - H2_" << j << "+=: " << H2_j << std::endl;
			#endif
									} // done sum over all quadrature points
									// if small, should be zero
									if ( fabs( H2_j ) < eps ) {
										H2_j = 0.0;
									}
									// add to rhs
									localH2rhs[ j ] += H2_j;
			#ifndef NLOG
									H2output = false;
									debugStream.Suspend(); // disable logging
			#endif
								} // done computing H2's volume integral
							}

							//                                                // we will call this one
							// (E)_{i,j} += -\int_{T}v_{j}\cdot\nabla q_{i}dx // E's volume integral
							//                                                // see also "E's entitity surface integral", "E's neighbour surface integral" and "E's boundary integral" below
			#ifndef NLOG
							if ( Edebug ) {
								debugStream.Resume(); // enable logging
								debugStream << "    = E volume =================" << std::endl;
								debugStream.Suspend();
							}
			#endif
							for ( int i = 0; i < numPressureBaseFunctionsElement; ++i ) {
								for ( int j = 0; j < numVelocityBaseFunctionsElement; ++j ) {
									double E_i_j = 0.0;
			#ifndef NLOG
			//                        if ( ( i == logBaseI ) && ( j == logBaseJ ) ) Eoutput = true;
			//                        if ( allOutput ) Eoutput = true;
									if ( Edebug ) Eoutput = true;
									if ( entityOutput && Eoutput ) debugStream.Resume(); // enable logging
									debugStream << "    = E ========================" << std::endl;
									debugStream << "    basefunctions " << i << " " << j << std::endl;
									debugStream << "    volumeQuadrature.nop() " << volumeQuadratureElement.nop() << std::endl;
			#endif
									// sum over all quadratur points
									for ( size_t quad = 0; quad < volumeQuadratureElement.nop(); ++ quad ) {
										// get x
										ElementCoordinateType x = volumeQuadratureElement.point( quad );
										// get the integration factor
										double elementVolume = geometry.integrationElement( x );
										// get the quadrature weight
										double integrationWeight = volumeQuadratureElement.weight( quad );
										// compute v_{j}\cdot(\nabla q_i)
										VelocityRangeType v_j( 0.0 );
										velocityBaseFunctionSetElement.evaluate( j, x, v_j );
										PressureJacobianRangeType jacobian_of_q_i( 0.0 );
										pressureBaseFunctionSetElement.jacobian( i, x, jacobian_of_q_i );
										const VelocityRangeType gradient_of_q_i_untransposed( jacobian_of_q_i[0] );
										const JacobianInverseTransposedType jacobianInverseTransposed = geometry.jacobianInverseTransposed( x );
										VelocityRangeType gradient_of_q_i( 0.0 );
										jacobianInverseTransposed.mv( gradient_of_q_i_untransposed, gradient_of_q_i );
										const double gradient_of_q_i_times_v_j = gradient_of_q_i * v_j;
										E_i_j += -1.0
											* elementVolume
											* integrationWeight
											* gradient_of_q_i_times_v_j;
			#ifndef NLOG
										debugStream << "    - quadPoint " << quad;
										Stuff::printFieldVector( x, "x", debugStream, "      " );
										debugStream << "\n        - elementVolume: " << elementVolume << std::endl;
										debugStream << "        - integrationWeight: " << integrationWeight;
										Stuff::printFieldVector( gradient_of_q_i, "gradient_of_q_i", debugStream, "      " );
										Stuff::printFieldVector( v_j, "v_j", debugStream, "      " );
										debugStream << "\n        - gradient_of_q_i_times_v_j: " << gradient_of_q_i_times_v_j << std::endl;
										debugStream << "        - E_" << i << "_" << j << "+=: " << E_i_j << std::endl;
			#endif
									} // done sum over all quadrature points
									// if small, should be zero
									if ( fabs( E_i_j ) < eps ) {
										E_i_j = 0.0;
									}
									else
										// add to matrix
										localEmatrixElement.add( i, j, E_i_j );
			#ifndef NLOG
									Eoutput = false;
									debugStream.Suspend(); // disable logging
			#endif
								}
							} // done computing E's volume integral

							// walk the intersections
							IntersectionIteratorType intItEnd = gridPart_.iend( entity );
							for (   IntersectionIteratorType intIt = gridPart_.ibegin( entity );
									intIt != intItEnd;
									++intIt ) {
										const typename IntersectionIteratorType::Intersection& intersection = *intIt;
			#ifndef NLOG
			//                    if ( ( outputIntersection == intersectionNR ) && entityOutput ) intersectionOutput = true;
								if ( entityOutput ) intersectionOutput = true;
								if ( intersectionOutput ) debugStream.Resume(); // enable logging
								debugStream << "    - ==== start calculations on intersection " << intersectionNR << std::endl;
			#endif

								// get intersection geometry
								typedef typename IntersectionIteratorType::Geometry
									IntersectionGeometryType;
								const IntersectionGeometryType& intersectionGeoemtry = intersection.intersectionGlobal();
			#ifndef NLOG
								// get corners
								debugStream << "      - intersection " << intersectionNR << " has " << intersectionGeoemtry.corners() << " corners:";
								for ( int i = 0; i < intersectionGeoemtry.corners(); ++i ) {
									const VelocityRangeType corner = intersectionGeoemtry.corner( i );
									Stuff::printFieldVector( corner, "corner_"+Stuff::toString(i), debugStream, "        " );
								}
								debugStream << "\n        length of intersection " << intersectionNR << " is " << Stuff::getLenghtOfIntersection( intersection ) << std::endl;
								debugStream.Suspend(); // disable logging
			#endif

								// get intersection quadrature, seen from inside
								const FaceQuadratureType faceQuadratureElement( gridPart_,
																				intersection,
																				( 4 * pressureSpaceOrder ) + 1,
																				FaceQuadratureType::INSIDE );

								// get flux coefficients
								const double lengthOfIntersection = Stuff::getLenghtOfIntersection( intersection );
								StabilizationCoefficients stabil_coeff ( discreteModel_.getStabilizationCoefficients() );
								const double C_11 = stabil_coeff.Factor("C11") * std::pow( lengthOfIntersection, stabil_coeff.Power("C11") );
								const double D_11 = stabil_coeff.Factor("D11") * std::pow( lengthOfIntersection, stabil_coeff.Power("D11") );
								//we'll leave this on 0 for the time being so it does not generate any additional  penalty terms
			//					const VelocityRangeType D_12(stabil_coeff.Factor("D12") );//TODO FIXME
								const VelocityRangeType D_12( 0 );//TODO FIXME

								// if we are inside the grid
								if ( intersection.neighbor() && !intersection.boundary() ) {
									// get neighbour
									const typename IntersectionIteratorType::EntityPointer neighbourPtr = intersection.outside();
									const EntityType& neighbour = *neighbourPtr;

									// get local matrices for the surface integrals
									LocalMInversMatrixType localMInversMatrixNeighbour = MInversMatrix.localMatrix( entity, neighbour );
									LocalWmatrixType localWmatrixNeighbour = Wmatrix.localMatrix( neighbour, entity );
									LocalXmatrixType localXmatrixNeighbour = Xmatrix.localMatrix( entity, neighbour );
									LocalYmatrixType localYmatrixNeighbour = Ymatrix.localMatrix( neighbour, entity );
									LocalZmatrixType localZmatrixNeighbour = Zmatrix.localMatrix( entity, neighbour );
									LocalEmatrixType localEmatrixNeighbour = Ematrix.localMatrix( neighbour, entity );
									LocalRmatrixType localRmatrixNeighbour = Rmatrix.localMatrix( entity, neighbour );

									// get basefunctionsets
									const SigmaBaseFunctionSetType sigmaBaseFunctionSetNeighbour = sigmaSpace_.baseFunctionSet( neighbour );
									const VelocityBaseFunctionSetType velocityBaseFunctionSetNeighbour = velocitySpace_.baseFunctionSet( neighbour );
									const PressureBaseFunctionSetType pressureBaseFunctionSetNeighbour = pressureSpace_.baseFunctionSet( neighbour );
									const int numSigmaBaseFunctionsNeighbour = sigmaBaseFunctionSetNeighbour.numBaseFunctions();
									const int numVelocityBaseFunctionsNeighbour = velocityBaseFunctionSetNeighbour.numBaseFunctions();
									const int numPressureBaseFunctionsNeighbour = pressureBaseFunctionSetNeighbour.numBaseFunctions();

									// get intersection quadrature, seen from outside
									const FaceQuadratureType faceQuadratureNeighbour(   gridPart_,
																						intersection,
																						( 4 * pressureSpaceOrder ) + 1,
																						FaceQuadratureType::OUTSIDE );

									// compute surface integrals

									//                                                                                                               // we will call this one
									// (W)_{i,j} += \int_{\varepsilon\in \Epsilon_{I}^{T}}-\hat{u}_{\sigma}^{U^{+}}(v_{j})\cdot\tau_{i}\cdot n_{T}ds // W's element surface integral
									//           += \int_{\varepsilon\in \Epsilon_{I}^{T}}-\hat{u}_{\sigma}^{U^{-}}(v_{j})\cdot\tau_{i}\cdot n_{T}ds // W's neighbour surface integral
									//                                                                                                               // see also "W's boundary integral" below
									//                                                                                                               // and "W's volume integral" above
			#ifndef NLOG
									if ( Wdebug ) {
										debugStream.Resume(); // enable logging
										debugStream << "      = W surface =======================" << std::endl;
										debugStream.Suspend();
									}
			#endif
			//                        if ( discreteModel_.hasVelocitySigmaFlux() ) {
										for ( int j = 0; j < numVelocityBaseFunctionsElement; ++j ) {
											// compute W's element surface integral
											for ( int i = 0; i < numSigmaBaseFunctionsElement; ++i ) {
												double W_i_j = 0.0;
			#ifndef NLOG
			//                                    if ( ( i == logBaseI ) && ( j == logBaseJ ) ) Woutput = true;
			//                                    if ( ( ( i == logBaseI ) && ( j == logBaseJ ) ) && Wdebug ) Woutput = true;
			//                                    if ( allOutput ) Woutput = true;
												if ( Wdebug ) Woutput = true;
												if ( intersectionOutput && Woutput ) debugStream.Resume(); // enable logging
												debugStream << "      = W element ======================" << std::endl;
												debugStream << "      basefunctions " << i << " " << j << std::endl;
												debugStream << "      faceQuadratureElement.nop() " << faceQuadratureElement.nop() << std::endl;
			#endif
												// sum over all quadrature points
												for ( size_t quad = 0; quad < faceQuadratureElement.nop(); ++quad ) {
													// get x in codim<0> and codim<1> coordinates
													const ElementCoordinateType x = faceQuadratureElement.point( quad );
													const LocalIntersectionCoordinateType xLocal = faceQuadratureElement.localPoint( quad );
													// get the integration factor
													const double elementVolume = intersectionGeoemtry.integrationElement( xLocal );
													// get the quadrature weight
													const double integrationWeight = faceQuadratureElement.weight( quad );
													// compute \hat{u}_{\sigma}^{U^{+}}(v_{j})\cdot\tau_{j}\cdot n_{T}
													const VelocityRangeType outerNormal = intersection.unitOuterNormal( xLocal );
													SigmaRangeType tau_i( 0.0 );
													sigmaBaseFunctionSetElement.evaluate( i, x, tau_i );
													VelocityRangeType v_j( 0.0 );
													velocityBaseFunctionSetElement.evaluate( j, x, v_j );
													VelocityRangeType tau_i_times_normal( 0.0 );
													tau_i.mv( outerNormal, tau_i_times_normal );
													double v_j_times_tau_i_times_normal = v_j * tau_i_times_normal;
													W_i_j += -0.5
														* elementVolume
														* integrationWeight
														* v_j_times_tau_i_times_normal;
			#ifndef NLOG
													debugStream << "      - quadPoint " << quad;
													Stuff::printFieldVector( x, "x", debugStream, "        " );
													Stuff::printFieldVector( xLocal, "xLocal", debugStream, "          " );
													Stuff::printFieldVector( outerNormal, "outerNormal", debugStream, "          " );
													debugStream << "\n        - elementVolume: " << elementVolume << std::endl;
													debugStream << "        - integrationWeight: " << integrationWeight;
													Stuff::printFieldMatrix( tau_i, "tau_i", debugStream, "        " );
													Stuff::printFieldVector( v_j, "v_j", debugStream, "        " );
													Stuff::printFieldVector( tau_i_times_normal, "tau_i_times_normal", debugStream, "        " );
													debugStream << "\n          - v_j_times_tau_i_times_normal: " << v_j_times_tau_i_times_normal << std::endl;
													debugStream << "          - W_" << i << "_" << j << "+=: " << W_i_j << std::endl;
			#endif
												} // done sum over all quadrature points
												// if small, should be zero
												if ( fabs( W_i_j ) < eps ) {
													W_i_j = 0.0;
												}
												else
													// add to matrix
													localWmatrixElement.add( i, j, W_i_j );
			#ifndef NLOG
												Woutput = false;
												debugStream.Suspend(); // disable logging
			#endif
											} // done computing W's element surface integral
											// compute W's neighbour surface integral
											for ( int i = 0; i < numSigmaBaseFunctionsNeighbour; ++i ) {
												double W_i_j = 0.0;
			#ifndef NLOG
			//                                    if ( ( i == logBaseI ) && ( j == logBaseJ ) ) Woutput = true;
			//                                    if ( ( ( i == logBaseI ) && ( j == logBaseJ ) ) && Wdebug ) Woutput = true;
			//                                    if ( allOutput ) Woutput = true;
												if ( Wdebug ) Woutput = true;
												if ( intersectionOutput && Woutput ) debugStream.Resume(); // enable logging
												debugStream << "      = W neighbour ====================" << std::endl;
												debugStream << "      basefunctions " << i << " " << j << std::endl;
												debugStream << "      faceQuadratureNeighbour.nop() " << faceQuadratureNeighbour.nop() << std::endl;
			#endif
												// sum over all quadrature points
												for ( size_t quad = 0; quad < faceQuadratureElement.nop(); ++quad ) {
													// get x in codim<0> and codim<1> coordinates
													const ElementCoordinateType xInside = faceQuadratureElement.point( quad );
													const ElementCoordinateType xOutside = faceQuadratureNeighbour.point( quad );
													const LocalIntersectionCoordinateType xLocal = faceQuadratureElement.localPoint( quad );
													// get the integration factor
													const double elementVolume = intersectionGeoemtry.integrationElement( xLocal );
													// get the quadrature weight
													const double integrationWeight = faceQuadratureElement.weight( quad );
													// compute \hat{u}_{\sigma}^{U^{-}}(v_{j})\cdot\tau_{j}\cdot n_{T}
													const VelocityRangeType outerNormal = intersection.unitOuterNormal( xLocal );
													SigmaRangeType tau_i( 0.0 );
													sigmaBaseFunctionSetNeighbour.evaluate( i, xOutside, tau_i );
													VelocityRangeType v_j( 0.0 );
													velocityBaseFunctionSetElement.evaluate( j, xInside, v_j );
													VelocityRangeType tau_i_times_normal( 0.0 );
													tau_i.mv( outerNormal, tau_i_times_normal );
													double v_j_times_tau_i_times_normal = v_j * tau_i_times_normal;
													W_i_j += 0.5
														* elementVolume
														* integrationWeight
														* v_j_times_tau_i_times_normal;
			#ifndef NLOG
													debugStream << "      - quadPoint " << quad;
													Stuff::printFieldVector( xInside, "xInside", debugStream, "        " );
													Stuff::printFieldVector( xOutside, "xOutside", debugStream, "        " );
													Stuff::printFieldVector( xLocal, "xLocal", debugStream, "          " );
													Stuff::printFieldVector( outerNormal, "outerNormal", debugStream, "          " );
													debugStream << "\n        - elementVolume: " << elementVolume << std::endl;
													debugStream << "        - integrationWeight: " << integrationWeight;
													Stuff::printFieldMatrix( tau_i, "tau_i", debugStream, "        " );
													Stuff::printFieldVector( v_j, "v_j", debugStream, "        " );
													Stuff::printFieldVector( tau_i_times_normal, "tau_i_times_normal", debugStream, "        " );
													debugStream << "\n          - v_j_times_tau_i_times_normal: " << v_j_times_tau_i_times_normal << std::endl;
													debugStream << "          - W_" << i << "_" << j << "+=: " << W_i_j << std::endl;
			#endif
												} // done sum over all quadrature points
												// if small, should be zero
												if ( fabs( W_i_j ) < eps ) {
													W_i_j = 0.0;
												}
												else
													// add to matrix
													localWmatrixNeighbour.add( i, j, W_i_j );
			#ifndef NLOG
												Woutput = false;
												debugStream.Suspend(); // disable logging
			#endif
											} // done computing W's neighbour surface integral
										} // done computing W's surface integrals
			//                        }


									//                                                                                                                   // we will call this one
									// (X)_{i,j} += \int_{\varepsilon\in\Epsilon_{I}^{T}}-\mu v_{i}\cdot\hat{\sigma}^{\sigma^{+}}(\tau_{j})\cdot n_{t}ds // X's element sourface integral
									//           += \int_{\varepsilon\in\Epsilon_{I}^{T}}-\mu v_{i}\cdot\hat{\sigma}^{\sigma^{-}}(\tau_{j})\cdot n_{t}ds // X's neighbour sourface integral
									//                                                                                                                   // see also "X's boundary integral" below
									//                                                                                                                   // and "X's volume integral" above
			#ifndef NLOG
									if ( Xdebug ) {
										debugStream.Resume(); // enable logging
										debugStream << "      = X surface =======================" << std::endl;
										debugStream.Suspend();
									}
			#endif
			//                        if ( discreteModel_.hasSigmaFlux() ) {
										for ( int j = 0; j < numSigmaBaseFunctionsElement; ++j ) {
											// compute X's element sourface integral
											for ( int i = 0; i < numVelocityBaseFunctionsElement; ++i ) {
												double X_i_j = 0.0;
			#ifndef NLOG
			//                                    if ( ( i == logBaseI ) && ( j == logBaseJ ) ) Xoutput = true;
			//                                    if ( allOutput ) Xoutput = true;
												if ( Xdebug ) Xoutput = true;
												if ( intersectionOutput && Xoutput ) debugStream.Resume(); // enable logging
												debugStream << "      = X element ======================" << std::endl;
												debugStream << "      basefunctions " << i << " " << j << std::endl;
												debugStream << "      volumeQuadratureElement.nop() " << volumeQuadratureElement.nop() << std::endl;
			#endif
												// sum over all quadrature points
												for ( size_t quad = 0; quad < faceQuadratureElement.nop(); ++quad ) {
													// get x in codim<0> and codim<1> coordinates
													const ElementCoordinateType x = faceQuadratureElement.point( quad );
													const LocalIntersectionCoordinateType xLocal = faceQuadratureElement.localPoint( quad );
													// get the integration factor
													const double elementVolume = intersectionGeoemtry.integrationElement( xLocal );
													// get the quadrature weight
													const double integrationWeight = faceQuadratureElement.weight( quad );
													// compute -\mu v_{i}\cdot\hat{\sigma}^{\sigma^{+}}(\tau_{j})\cdot n_{t}
													const VelocityRangeType outerNormal = intersection.unitOuterNormal( xLocal );
													SigmaRangeType tau_j( 0.0 );
													sigmaBaseFunctionSetElement.evaluate( j, x, tau_j );
													VelocityRangeType tau_j_times_normal( 0.0 );
													tau_j.mv( outerNormal, tau_j_times_normal );
													const double v_i_times_tau_j_times_normal = velocityBaseFunctionSetElement.evaluateSingle( i, x, tau_j_times_normal );
													X_i_j += -0.5
														* elementVolume
														* integrationWeight
														* mu
														* v_i_times_tau_j_times_normal;
			#ifndef NLOG
													debugStream << "      - quadPoint " << quad;
													Stuff::printFieldVector( x, "x", debugStream, "        " );
													Stuff::printFieldVector( xLocal, "xLocal", debugStream, "          " );
													Stuff::printFieldVector( outerNormal, "outerNormal", debugStream, "          " );
													debugStream << "\n        - elementVolume: " << elementVolume << std::endl;
													debugStream << "        - integrationWeight: " << integrationWeight;
													Stuff::printFieldMatrix( tau_j, "tau_j", debugStream, "        " );
													Stuff::printFieldVector( tau_j_times_normal, "tau_j_times_normal", debugStream, "        " );
													debugStream << "\n          - v_i_times_tau_j_times_normal: " << v_i_times_tau_j_times_normal << std::endl;
													debugStream << "          - X_" << i << "_" << j << "+=: " << X_i_j << std::endl;
			#endif
												} // done sum over all quadrature points
												// if small, should be zero
												if ( fabs( X_i_j ) < eps ) {
													X_i_j = 0.0;
												}
												else
													// add to matrix
													localXmatrixElement.add( i, j, X_i_j );
			#ifndef NLOG
												Xoutput = false;
												debugStream.Suspend(); // disable logging
			#endif
											} // done computing X's element sourface integral
											// compute X's neighbour sourface integral
											for ( int i = 0; i < numVelocityBaseFunctionsNeighbour; ++i ) {
												double X_i_j = 0.0;
			#ifndef NLOG
			//                                    if ( ( i == logBaseI ) && ( j == logBaseJ ) ) Xoutput = true;
			//                                    if ( allOutput ) Xoutput = true;
												if ( Xdebug ) Xoutput = true;
												if ( intersectionOutput && Xoutput ) debugStream.Resume(); // enable logging
												debugStream << "      = X neighbour ====================" << std::endl;
												debugStream << "      basefunctions " << i << " " << j << std::endl;
												debugStream << "      faceQuadratureNeighbour.nop() " << faceQuadratureNeighbour.nop() << std::endl;
			#endif
												// sum over all quadrature points
												for ( size_t quad = 0; quad < faceQuadratureNeighbour.nop(); ++quad ) {
													// get x codim<0> and codim<1> coordinates
													const ElementCoordinateType xInside = faceQuadratureElement.point( quad );
													const ElementCoordinateType xOutside = faceQuadratureNeighbour.point( quad );
													const LocalIntersectionCoordinateType xLocal = faceQuadratureNeighbour.localPoint( quad );
													// get the integration factor
													const double elementVolume = intersectionGeoemtry.integrationElement( xLocal );
													// get the quadrature weight
													const double integrationWeight = faceQuadratureElement.weight( quad );
													// compute -\mu v_{i}\cdot\hat{\sigma}^{\sigma^{-}}(\tau_{j})\cdot n_{t}
													const VelocityRangeType outerNormal = intersection.unitOuterNormal( xLocal );
													SigmaRangeType tau_j( 0.0 );
													sigmaBaseFunctionSetNeighbour.evaluate( j, xOutside, tau_j );
													VelocityRangeType tau_j_times_normal( 0.0 );
													tau_j.mv( outerNormal, tau_j_times_normal );
													const double v_i_times_tau_j_times_normal = velocityBaseFunctionSetElement.evaluateSingle( i, xInside, tau_j_times_normal );
													X_i_j += -0.5
														* elementVolume
														* integrationWeight
														* mu
														* v_i_times_tau_j_times_normal;
			#ifndef NLOG
													debugStream << "      - quadPoint " << quad;
													Stuff::printFieldVector( xInside, "xInside", debugStream, "        " );
													Stuff::printFieldVector( xOutside, "xOutside", debugStream, "        " );
													Stuff::printFieldVector( xLocal, "xLocal", debugStream, "          " );
													Stuff::printFieldVector( outerNormal, "outerNormal", debugStream, "          " );
													debugStream << "\n        - elementVolume: " << elementVolume << std::endl;
													debugStream << "        - integrationWeight: " << integrationWeight;
													Stuff::printFieldMatrix( tau_j, "tau_j", debugStream, "        " );
													Stuff::printFieldVector( tau_j_times_normal, "tau_j_times_normal", debugStream, "        " );
													debugStream << "\n          - v_i_times_tau_j_times_normal: " << v_i_times_tau_j_times_normal << std::endl;
													debugStream << "          - X_" << i << "_" << j << "+=: " << X_i_j << std::endl;
			#endif
												} // done sum over all quadrature points
												// if small, should be zero
												if ( fabs( X_i_j ) < eps ) {
													X_i_j = 0.0;
												}
												else
													// add to matrix
													localXmatrixNeighbour.add( i, j, X_i_j );
			#ifndef NLOG
												Xoutput = false;
												debugStream.Suspend(); // disable logging
			#endif
											} // done computing X's neighbour sourface integral
										} // done computing X's sourface integrals
			//                        }

									//                                                                                                         // we call this one
									// (Y)_{i,j} += \int_{\varepsilon\in\Epsilon_{I}^{T}}-\mu v_{i}\cdot\hat{\sigma}^{U{+}}(v{j})\cdot n_{t}ds // Y's element surface integral
									//           += \int_{\varepsilon\in\Epsilon_{I}^{T}}-\mu v_{i}\cdot\hat{\sigma}^{U{-}}(v{j})\cdot n_{t}ds // Y's neighbour surface integral
									//                                                                                                         // see also "Y's boundary integral" below
			#ifndef NLOG
									if ( Ydebug ) {
										debugStream.Resume(); // enable logging
										debugStream << "      = Y surface =======================" << std::endl;
										debugStream.Suspend();
									}
			#endif
			//                        if ( discreteModel_.hasSigmaFlux() ) {
										for ( int j = 0; j < numVelocityBaseFunctionsElement; ++j ) {
											// compute Y's element surface integral
											for ( int i = 0; i < numVelocityBaseFunctionsElement; ++i ) {
												double Y_i_j = 0.0;
			#ifndef NLOG
			//                                    if ( ( i == logBaseI ) && ( j == logBaseJ ) ) Youtput = true;
			//                                    if ( allOutput ) Youtput = true;
												if ( Ydebug ) Youtput = true;
												if ( intersectionOutput && Youtput ) debugStream.Resume(); // enable logging
												debugStream << "      = Y element ======================" << std::endl;
												debugStream << "      basefunctions " << i << " " << j << std::endl;
												debugStream << "      faceQuadratureElement.nop() " << faceQuadratureElement.nop() << std::endl;
			#endif
												// sum over all quadrature points
												for ( size_t quad = 0; quad < faceQuadratureElement.nop(); ++quad ) {
													// get x codim<0> and codim<1> coordinates
													const ElementCoordinateType x = faceQuadratureElement.point( quad );
													const ElementCoordinateType xInside = faceQuadratureElement.point( quad );
													const ElementCoordinateType xOutside = faceQuadratureNeighbour.point( quad );
													const LocalIntersectionCoordinateType xLocal = faceQuadratureElement.localPoint( quad );
													const VelocityRangeType xWorld = geometry.global( x );
													const VelocityRangeType xWorld_Outside = geometry.global( xOutside );
													// get the integration factor
													const double elementVolume = intersectionGeoemtry.integrationElement( xLocal );
													// get the quadrature weight
													const double integrationWeight = faceQuadratureElement.weight( quad );
													// compute -\mu v_{i}\cdot\hat{\sigma}^{U{+}}(v{j})\cdot n_{t}
													const VelocityRangeType outerNormal = intersection.unitOuterNormal( xLocal );
													VelocityRangeType v_j( 0.0 );
													velocityBaseFunctionSetElement.evaluate( j, x, v_j );
													VelocityRangeType v_i( 0.0 );
													velocityBaseFunctionSetElement.evaluate( i, x, v_i );
													const double v_i_times_v_j = v_i * v_j;
													Y_i_j += C_11
														* elementVolume
														* integrationWeight
														* mu
														* v_i_times_v_j;
													//calc u^c_h \tensor beta * v \tensor n
													VelocityRangeType beta_eval;
													beta_.evaluate( xWorld, beta_eval );
													const double beta_times_normal = beta_eval * outerNormal;
													VelocityJacobianRangeType v_i_tensor_n = dyadicProduct( v_i, outerNormal );
													double c_s;

													if ( beta_times_normal < 0 ) {
														c_s = - beta_times_normal * 0.5;
													}
													else {
														c_s = beta_times_normal * 0.5;
													}
													VelocityRangeType u_h_half = v_i;
													u_h_half *= 0.5;
													VelocityJacobianRangeType flux_value = dyadicProduct( v_i, beta_eval );

													VelocityJacobianRangeType u_jump = dyadicProduct( v_i, outerNormal );
													u_jump *= c_s;

													flux_value += u_jump;

													double ret  = Stuff::colonProduct( flux_value, v_i_tensor_n );
													//inner edge (self)
													Y_i_j += elementVolume
														* integrationWeight
														* ret;

			#ifndef NLOG
													debugStream << "      - quadPoint " << quad;
													Stuff::printFieldVector( x, "x", debugStream, "        " );
													Stuff::printFieldVector( xLocal, "xLocal", debugStream, "          " );
													Stuff::printFieldVector( outerNormal, "outerNormal", debugStream, "          " );
													debugStream << "\n        - elementVolume: " << elementVolume << std::endl;
													debugStream << "        - integrationWeight: " << integrationWeight;
													Stuff::printFieldVector( v_i, "v_i", debugStream, "        " );
													Stuff::printFieldVector( v_j, "v_j", debugStream, "        " );
													debugStream << "\n          - v_i_times_v_j: " << v_i_times_v_j << std::endl;
													debugStream << "          - Y_" << i << "_" << j << "+=: " << Y_i_j << std::endl;
			#endif
												} // done sum over all quadrature points
												// if small, should be zero
												if ( fabs( Y_i_j ) < eps ) {
													Y_i_j = 0.0;
												}
												else
													// add to matrix
													localYmatrixElement.add( i, j, Y_i_j );
			#ifndef NLOG
												Youtput = false;
												debugStream.Suspend(); // disable logging
			#endif
											} // done computing Y's element surface integral
											// compute Y's neighbour surface integral
											for ( int i = 0; i < numVelocityBaseFunctionsNeighbour; ++i ) {
												double Y_i_j = 0.0;
			#ifndef NLOG
			//                                    if ( ( i == logBaseI ) && ( j == logBaseJ ) ) Youtput = true;
			//                                    if ( allOutput ) Youtput = true;
												if ( Ydebug ) Youtput = true;
												if ( intersectionOutput && Youtput ) debugStream.Resume(); // enable logging
												debugStream << "      = Y neighbour ====================" << std::endl;
												debugStream << "      basefunctions " << i << " " << j << std::endl;
												debugStream << "      faceQuadratureNeighbour.nop() " << faceQuadratureNeighbour.nop() << std::endl;
			#endif
												// sum over all quadrature points
												for ( size_t quad = 0; quad < faceQuadratureNeighbour.nop(); ++quad ) {
													// get x codim<0> and codim<1> coordinates
													const ElementCoordinateType xInside = faceQuadratureElement.point( quad );
													const ElementCoordinateType xOutside = faceQuadratureNeighbour.point( quad );
													const VelocityRangeType xWorld = geometry.global( xInside );
													const VelocityRangeType xWorld_Outside = geometry.global( xOutside );
													const LocalIntersectionCoordinateType xLocal = faceQuadratureNeighbour.localPoint( quad );
													// get the integration factor
													const double elementVolume = intersectionGeoemtry.integrationElement( xLocal );
													// get the quadrature weight
													const double integrationWeight = faceQuadratureNeighbour.weight( quad );
													// compute -\mu v_{i}\cdot\hat{\sigma}^{U{-}}(v{j})\cdot n_{t}
													const VelocityRangeType outerNormal = intersection.unitOuterNormal( xLocal );
													VelocityRangeType v_i( 0.0 );
													velocityBaseFunctionSetNeighbour.evaluate( i, xOutside, v_i );
													VelocityRangeType v_j( 0.0 );
													velocityBaseFunctionSetElement.evaluate( j, xInside, v_j );
													const double v_i_times_v_j = v_i * v_j;
													Y_i_j += -1.0
														* C_11
														* elementVolume
														* integrationWeight
														* mu
														* v_i_times_v_j;
													//calc u^c_h \tensor beta * v \tensor n
													VelocityRangeType beta_eval;
													beta_.evaluate( xWorld, beta_eval );
													const double beta_times_normal = beta_eval * outerNormal;
													VelocityJacobianRangeType v_i_tensor_n = dyadicProduct( v_i, outerNormal );
													double c_s;

													if ( beta_times_normal < 0 ) {
														c_s = - beta_times_normal * 0.5;
													}
													else {
														c_s = beta_times_normal * 0.5;
													}
													VelocityRangeType u_h_half = v_i;
													u_h_half *= 0.5;
													VelocityJacobianRangeType flux_value = dyadicProduct( v_i, beta_eval );

													VelocityJacobianRangeType u_jump = dyadicProduct( v_i, outerNormal );
													u_jump *= c_s;

													flux_value += u_jump;

													double ret  = Stuff::colonProduct( flux_value, v_i_tensor_n );
													//inner edge (self)
													Y_i_j += elementVolume
														* integrationWeight
														* ret;
			#ifndef NLOG
													debugStream << "      - quadPoint " << quad;
													Stuff::printFieldVector( xInside, "xInside", debugStream, "        " );
													Stuff::printFieldVector( xOutside, "xOutside", debugStream, "        " );
													Stuff::printFieldVector( xLocal, "xLocal", debugStream, "          " );
													Stuff::printFieldVector( outerNormal, "outerNormal", debugStream, "          " );
													debugStream << "\n        - elementVolume: " << elementVolume << std::endl;
													debugStream << "        - integrationWeight: " << integrationWeight;
													Stuff::printFieldVector( v_i, "v_i", debugStream, "        " );
													Stuff::printFieldVector( v_j, "v_j", debugStream, "        " );
													debugStream << "\n          - v_i_times_v_j: " << v_i_times_v_j << std::endl;
													debugStream << "          - Y_" << i << "_" << j << "+=: " << Y_i_j << std::endl;
			#endif
												} // done sum over all quadrature points
												// if small, should be zero
												if ( fabs( Y_i_j ) < eps ) {
													Y_i_j = 0.0;
												}
												else
													// add to matrix
													localYmatrixNeighbour.add( i, j, Y_i_j );
			#ifndef NLOG
												Youtput = false;
												debugStream.Suspend(); // disable logging
			#endif
											} // done computing Y's neighbour surface integral
										} // done computing Y's surface integrals
			//                        }

									//                                                                                                  // we will call this one
									// (Z)_{i,j} += \int_{\varepsilon\in\Epsilon_{I}^{T}}\hat{p}^{P^{+}}(q_{j})\cdot v_{i}\cdot n_{T}ds // Z's element surface integral
									//           += \int_{\varepsilon\in\Epsilon_{I}^{T}}\hat{p}^{P^{-}}(q_{j})\cdot v_{i}\cdot n_{T}ds // Z's neighbour surface integral
									//                                                                                                  // see also "Z's boundary integral" below
									//                                                                                                  // and "Z's volume integral" above
			#ifndef NLOG
									if ( Zdebug ) {
										debugStream.Resume(); // enable logging
										debugStream << "      = Z surface =======================" << std::endl;
										debugStream.Suspend();
									}
			#endif
			//                        if ( discreteModel_.hasPressureFlux() ) {
										for ( int j = 0; j < numPressureBaseFunctionsElement; ++j ) {
											// compute Z's element surface integral
											for ( int i = 0; i < numVelocityBaseFunctionsElement; ++i ) {
												double Z_i_j = 0.0;
			#ifndef NLOG
			//                                    if ( ( i == logBaseI ) && ( j == logBaseJ ) ) Zoutput = true;
			//                                    if ( allOutput ) Zoutput = true;
												if ( Zdebug ) Zoutput = true;
												if ( intersectionOutput && Zoutput ) debugStream.Resume(); // enable logging
												debugStream << "      = Z element ======================" << std::endl;
												debugStream << "      basefunctions " << i << " " << j << std::endl;
												debugStream << "      faceQuadratureElement.nop() " << faceQuadratureElement.nop() << std::endl;
			#endif
												// sum over all quadrature points
												for ( size_t quad = 0; quad < faceQuadratureElement.nop(); ++quad ) {
													// get x codim<0> and codim<1> coordinates
													const ElementCoordinateType x = faceQuadratureElement.point( quad );
													const LocalIntersectionCoordinateType xLocal = faceQuadratureElement.localPoint( quad );
													// get the integration factor
													const double elementVolume = intersectionGeoemtry.integrationElement( xLocal );
													// get the quadrature weight
													const double integrationWeight = faceQuadratureElement.weight( quad );
													// compute \hat{p}^{P^{+}}(q_{j})\cdot v_{i}\cdot n_{T}
													const VelocityRangeType outerNormal = intersection.unitOuterNormal( xLocal );
													VelocityRangeType v_i( 0.0 );
													velocityBaseFunctionSetElement.evaluate( i, x, v_i );
													PressureRangeType q_j( 0.0 );
													pressureBaseFunctionSetElement.evaluate( j, x, q_j );
													const double v_i_times_normal = v_i * outerNormal;
													const double p_factor = ( 0.5 - ( D_12 * outerNormal ) );// (0.5 p - p D_12 ) n ) <- p+
			//										const double p_factor = ( 0.5 - ( 1 ) );// (0.5 p - p D_12 ) n ) <- p+
													const double q_j_times_v_i_times_normal =  q_j * v_i_times_normal;
													Z_i_j += p_factor
														* elementVolume
														* integrationWeight
														* q_j_times_v_i_times_normal;
			#ifndef NLOG
													debugStream << "      - quadPoint " << quad;
													Stuff::printFieldVector( x, "x", debugStream, "        " );
													Stuff::printFieldVector( xLocal, "xLocal", debugStream, "          " );
													Stuff::printFieldVector( outerNormal, "outerNormal", debugStream, "          " );
													debugStream << "\n        - elementVolume: " << elementVolume << std::endl;
													debugStream << "        - integrationWeight: " << integrationWeight;
													Stuff::printFieldVector( v_i, "v_i", debugStream, "        " );
													Stuff::printFieldVector( q_j, "q_j", debugStream, "        " );
													debugStream << "\n          - v_i_times_normal: " << v_i_times_normal << std::endl;
													debugStream << "          - q_j_times_v_i_times_normal: " << q_j_times_v_i_times_normal << std::endl;
													debugStream << "          - Z_" << i << "_" << j << "+=: " << Z_i_j << std::endl;
			#endif
												} // done sum over all quadrature points
												// if small, should be zero
												if ( fabs( Z_i_j ) < eps ) {
													Z_i_j = 0.0;
												}
												else
													// add to matrix
													localZmatrixElement.add( i, j, Z_i_j );
			#ifndef NLOG
												Zoutput = false;
												debugStream.Suspend(); // disable logging
			#endif
											} // done computing Z's element surface integral
											// compute Z's neighbour surface integral
											for ( int i = 0; i < numVelocityBaseFunctionsNeighbour; ++i ) {
												double Z_i_j = 0.0;
			#ifndef NLOG
			//                                    if ( ( i == logBaseI ) && ( j == logBaseJ ) ) Zoutput = true;
			//                                    if ( allOutput ) Zoutput = true;
												if ( Zdebug ) Zoutput = true;
												if ( intersectionOutput && Zoutput ) debugStream.Resume(); // enable logging
												debugStream << "      = Z neighbour ====================" << std::endl;
												debugStream << "      basefunctions " << i << " " << j << std::endl;
												debugStream << "      faceQuadratureNeighbour.nop() " << faceQuadratureNeighbour.nop() << std::endl;
			#endif
												// sum over all quadrature points
												for ( size_t quad = 0; quad < faceQuadratureNeighbour.nop(); ++quad ) {
													// get x codim<0> and codim<1> coordinates
													const ElementCoordinateType xInside = faceQuadratureElement.point( quad );
													const ElementCoordinateType xOutside = faceQuadratureNeighbour.point( quad );
													const LocalIntersectionCoordinateType xLocal = faceQuadratureNeighbour.localPoint( quad );
													// get the integration factor
													const double elementVolume = intersectionGeoemtry.integrationElement( xLocal );
													// get the quadrature weight
													const double integrationWeight = faceQuadratureNeighbour.weight( quad );
													// compute \hat{p}^{P^{+}}(q_{j})\cdot v_{i}\cdot n_{T}
													const VelocityRangeType outerNormal = intersection.unitOuterNormal( xLocal );
													VelocityRangeType v_i( 0.0 );
													velocityBaseFunctionSetElement.evaluate( i, xInside, v_i );
													PressureRangeType q_j( 0.0 );
													pressureBaseFunctionSetNeighbour.evaluate( j, xOutside, q_j );
													const double v_i_times_normal = v_i * outerNormal;
													const double p_factor = ( 0.5 + ( D_12 * outerNormal ) );// (0.5 p + p D_12 ) n ) <- p-
			//										const double p_factor = ( 0.5 + ( 1 ) );// (0.5 p + p D_12 ) n ) <- p-
													const double q_j_times_v_i_times_normal = q_j * v_i_times_normal;
													Z_i_j += p_factor
														* elementVolume
														* integrationWeight
														* q_j_times_v_i_times_normal;
			#ifndef NLOG
													debugStream << "      - quadPoint " << quad;
													Stuff::printFieldVector( xInside, "xInside", debugStream, "        " );
													Stuff::printFieldVector( xOutside, "xOutside", debugStream, "        " );
													Stuff::printFieldVector( xLocal, "xLocal", debugStream, "          " );
													Stuff::printFieldVector( outerNormal, "outerNormal", debugStream, "          " );
													debugStream << "\n        - elementVolume: " << elementVolume << std::endl;
													debugStream << "        - integrationWeight: " << integrationWeight;
													Stuff::printFieldVector( v_i, "v_i", debugStream, "        " );
													Stuff::printFieldVector( q_j, "q_j", debugStream, "        " );
													debugStream << "\n          - v_i_times_normal: " << v_i_times_normal << std::endl;
													debugStream << "          - q_j_times_v_i_times_normal: " << q_j_times_v_i_times_normal << std::endl;
													debugStream << "          - Z_" << i << "_" << j << "+=: " << Z_i_j << std::endl;
			#endif
												} // done sum over all quadrature points
												// if small, should be zero
												if ( fabs( Z_i_j ) < eps ) {
													Z_i_j = 0.0;
												}
												else
													// add to matrix
													localZmatrixNeighbour.add( i, j, Z_i_j );
			#ifndef NLOG
												Zoutput = false;
												debugStream.Suspend(); // disable logging
			#endif
											} // done computing Z's neighbour surface integral
										} // done computing Z's surface integrals
			//                        }

									//                                                                                                // we will call this one
									// (E)_{i,j} += \int_{\varepsilon\in\Epsilon_{I}^{T}}\hat{u}_{p}^{U^{+}}(v_{j})\cdot n_{T}q_{i}ds // E's element surface integral
									//           += \int_{\varepsilon\in\Epsilon_{I}^{T}}\hat{u}_{p}^{U^{-}}(v_{j})\cdot n_{T}q_{i}ds // E's neighbour surface integral
									//                                                                                                // see also "E's boundary integral" below
									//                                                                                                // and "E's volume integral" above
			#ifndef NLOG
									if ( Edebug ) {
										debugStream.Resume(); // enable logging
										debugStream << "      = E surface =======================" << std::endl;
										debugStream.Suspend();
									}
			#endif
			//                        if ( discreteModel_.hasVelocityPressureFlux() ) {
										for ( int j = 0; j < numVelocityBaseFunctionsElement; ++j ) {
											// compute E's element surface integral
											for ( int i = 0; i < numPressureBaseFunctionsElement; ++i ) {
												double E_i_j = 0.0;
			#ifndef NLOG
			//                                    if ( ( i == logBaseI ) && ( j == logBaseJ ) ) Eoutput = true;
			//                                    if ( allOutput ) Eoutput = true;
												if ( Edebug ) Eoutput = true;
												if ( intersectionOutput && Eoutput ) debugStream.Resume(); // enable logging
												debugStream << "      = E element ======================" << std::endl;
												debugStream << "      basefunctions " << i << " " << j << std::endl;
												debugStream << "      faceQuadratureElement.nop() " << faceQuadratureElement.nop() << std::endl;
			#endif
												// sum over all quadrature points
												for ( size_t quad = 0; quad < faceQuadratureElement.nop(); ++quad ) {
													// get x codim<0> and codim<1> coordinates
													const ElementCoordinateType x = faceQuadratureElement.point( quad );
													const LocalIntersectionCoordinateType xLocal = faceQuadratureElement.localPoint( quad );
													// get the integration factor
													const double elementVolume = intersectionGeoemtry.integrationElement( xLocal );
													// get the quadrature weight
													const double integrationWeight = faceQuadratureElement.weight( quad );
													// compute \hat{u}_{p}^{U^{+}}(v_{j})\cdot n_{T}q_{i}
													const VelocityRangeType outerNormal = intersection.unitOuterNormal( xLocal );
													VelocityRangeType v_j( 0.0 );
													velocityBaseFunctionSetElement.evaluate( j, x, v_j );
													PressureRangeType q_i( 0.0 );
													pressureBaseFunctionSetElement.evaluate( i, x, q_i );
													const double v_j_times_normal = v_j * outerNormal;
													const double q_i_times_v_j_times_normal = q_i * v_j_times_normal;
													E_i_j += 0.5
														* elementVolume
														* integrationWeight
														* q_i_times_v_j_times_normal;
			#ifndef NLOG
													debugStream << "      - quadPoint " << quad;
													Stuff::printFieldVector( x, "x", debugStream, "        " );
													Stuff::printFieldVector( xLocal, "xLocal", debugStream, "          " );
													Stuff::printFieldVector( outerNormal, "outerNormal", debugStream, "          " );
													debugStream << "\n        - elementVolume: " << elementVolume << std::endl;
													debugStream << "        - integrationWeight: " << integrationWeight;
													Stuff::printFieldVector( q_i, "q_i", debugStream, "        " );
													Stuff::printFieldVector( v_j, "v_j", debugStream, "        " );
													debugStream << "\n          - v_j_times_normal: " << v_j_times_normal << std::endl;
													debugStream << "          - q_i_times_v_j_times_normal: " << q_i_times_v_j_times_normal << std::endl;
													debugStream << "          - E_" << i << "_" << j << "+=: " << E_i_j << std::endl;
			#endif
												} // done sum over all quadrature points
												// if small, should be zero
												if ( fabs( E_i_j ) < eps ) {
													E_i_j = 0.0;
												}
												else
													// add to matrix
													localEmatrixElement.add( i, j, E_i_j );
			#ifndef NLOG
												Eoutput = false;
												debugStream.Suspend(); // disable logging
			#endif
											} // done computing E's element surface integral
											// compute E's neighbour surface integral
											for ( int i = 0; i < numPressureBaseFunctionsNeighbour; ++i ) {
												double E_i_j = 0.0;
			#ifndef NLOG
			//                                    if ( ( i == logBaseI ) && ( j == logBaseJ ) ) Eoutput = true;
			//                                    if ( allOutput ) Eoutput = true;
												if ( Edebug ) Eoutput = true;
												if ( intersectionOutput && Eoutput ) debugStream.Resume(); // enable logging
												debugStream << "      = E neighbour ====================" << std::endl;
												debugStream << "      basefunctions " << i << " " << j << std::endl;
												debugStream << "      faceQuadratureNeighbour.nop() " << faceQuadratureNeighbour.nop() << std::endl;
			#endif
												// sum over all quadrature points
												for ( size_t quad = 0; quad < faceQuadratureNeighbour.nop(); ++quad ) {
													// get x codim<0> and codim<1> coordinates
													const ElementCoordinateType xInside = faceQuadratureElement.point( quad );
													const ElementCoordinateType xOutside = faceQuadratureNeighbour.point( quad );
													const LocalIntersectionCoordinateType xLocal = faceQuadratureNeighbour.localPoint( quad );
													// get the integration factor
													const double elementVolume = intersectionGeoemtry.integrationElement( xLocal );
													// get the quadrature weight
													const double integrationWeight = faceQuadratureNeighbour.weight( quad );
													// compute \hat{u}_{p}^{U^{-}}(v_{j})\cdot n_{T}q_{i}
													const VelocityRangeType outerNormal = intersection.unitOuterNormal( xLocal );
													VelocityRangeType v_j( 0.0 );
													velocityBaseFunctionSetElement.evaluate( j, xInside, v_j );
													PressureRangeType q_i( 0.0 );
													pressureBaseFunctionSetNeighbour.evaluate( i, xOutside, q_i );
													const double v_j_times_normal = v_j * outerNormal;
													const double q_i_times_v_j_times_normal = q_i * v_j_times_normal;
													E_i_j += -0.5
														* elementVolume
														* integrationWeight
														* q_i_times_v_j_times_normal;
			#ifndef NLOG
													debugStream << "      - quadPoint " << quad;
													Stuff::printFieldVector( xInside, "xInside", debugStream, "        " );
													Stuff::printFieldVector( xOutside, "xOutside", debugStream, "        " );
													Stuff::printFieldVector( xLocal, "xLocal", debugStream, "          " );
													Stuff::printFieldVector( outerNormal, "outerNormal", debugStream, "          " );
													debugStream << "\n        - elementVolume: " << elementVolume << std::endl;
													debugStream << "        - integrationWeight: " << integrationWeight;
													Stuff::printFieldVector( q_i, "q_i", debugStream, "        " );
													Stuff::printFieldVector( v_j, "v_j", debugStream, "        " );
													debugStream << "\n          - v_j_times_normal: " << v_j_times_normal << std::endl;
													debugStream << "\n          - q_i_times_v_j_times_normal: " << q_i_times_v_j_times_normal << std::endl;
													debugStream << "          - E_" << i << "_" << j << "+=: " << E_i_j << std::endl;
			#endif
												} // done sum over all quadrature points
												// if small, should be zero
												if ( fabs( E_i_j ) < eps ) {
													E_i_j = 0.0;
												}
												else
													// add to matrix
													localEmatrixNeighbour.add( i, j, E_i_j );
			#ifndef NLOG
												Eoutput = false;
												debugStream.Suspend(); // disable logging
			#endif
											} // done computing E's neighbour surface integral
										} // done computing E's surface integrals
			//                        }

									//                                                                                                // we will call this one
									// (R)_{i,j} += \int_{\varepsilon\in\Epsilon_{I}^{T}}\hat{u}_{p}^{P^{+}}(q_{j})\cdot n_{T}q_{i}ds // R's element surface integral
									//           += \int_{\varepsilon\in\Epsilon_{I}^{T}}\hat{u}_{p}^{P^{-}}(q_{j})\cdot n_{T}q_{i}ds // R's neighbour surface integral
									//                                                                                                // see also "R's boundary integral" below
			#ifndef NLOG
									if ( Rdebug ) {
										debugStream.Resume(); // enable logging
										debugStream << "      = R surface =======================" << std::endl;
										debugStream.Suspend();
									}
			#endif
			//                        if ( discreteModel_.hasVelocityPressureFlux() ) {
										for ( int j = 0; j < numPressureBaseFunctionsElement; ++j ) {
											// compute R's element surface integral
											for ( int i = 0; i < numPressureBaseFunctionsElement; ++i ) {
												double R_i_j = 0.0;
			#ifndef NLOG
			//                                    if ( ( i == logBaseI ) && ( j == logBaseJ ) ) Routput = true;
			//                                    if ( allOutput ) Routput = true;
												if ( Rdebug ) Routput = true;
												if ( intersectionOutput && Routput ) debugStream.Resume(); // enable logging
												debugStream << "      = R element ======================" << std::endl;
												debugStream << "      basefunctions " << i << " " << j << std::endl;
												debugStream << "      faceQuadratureElement.nop() " << faceQuadratureElement.nop() << std::endl;
			#endif
												// sum over all quadrature points
												for ( size_t quad = 0; quad < faceQuadratureElement.nop(); ++quad ) {
													// get x codim<0> and codim<1> coordinates
													const ElementCoordinateType x = faceQuadratureElement.point( quad );
													const LocalIntersectionCoordinateType xLocal = faceQuadratureElement.localPoint( quad );
													// get the integration factor
													const double elementVolume = intersectionGeoemtry.integrationElement( xLocal );
													// get the quadrature weight
													const double integrationWeight = faceQuadratureElement.weight( quad );
													// compute \hat{u}_{p}^{P^{+}}(q_{j})\cdot n_{T}q_{i}
													const VelocityRangeType outerNormal = intersection.unitOuterNormal( xLocal );
													PressureRangeType q_i( 0.0 );
													pressureBaseFunctionSetElement.evaluate( i, x, q_i );
													PressureRangeType q_j( 0.0 );
													pressureBaseFunctionSetElement.evaluate( j, x, q_j );
													const double q_i_times_q_j = q_i * q_j;
													R_i_j += D_11
														* elementVolume
														* integrationWeight
														* q_i_times_q_j;
			#ifndef NLOG
													debugStream << "      - quadPoint " << quad;
													Stuff::printFieldVector( x, "x", debugStream, "        " );
													Stuff::printFieldVector( xLocal, "xLocal", debugStream, "          " );
													Stuff::printFieldVector( outerNormal, "outerNormal", debugStream, "          " );
													debugStream << "\n        - elementVolume: " << elementVolume << std::endl;
													debugStream << "        - integrationWeight: " << integrationWeight;
													Stuff::printFieldVector( q_i, "q_i", debugStream, "        " );
													Stuff::printFieldVector( q_j, "q_j", debugStream, "        " );
													debugStream << "\n          - q_i_times_q_j: " << q_i_times_q_j << std::endl;
			#endif
												} // done sum over all quadrature points
												// if small, should be zero
												if ( fabs( R_i_j ) < eps ) {
													R_i_j = 0.0;
												}
												else
													// add to matrix
													localRmatrixElement.add( i, j, R_i_j );
			#ifndef NLOG
												Routput = false;
												debugStream.Suspend(); // disable logging
			#endif
											} // done computing R's element surface integral
											// compute R's neighbour surface integral
											for ( int i = 0; i < numPressureBaseFunctionsNeighbour; ++i ) {
												double R_i_j = 0.0;
			#ifndef NLOG
			//                                    if ( ( i == logBaseI ) && ( j == logBaseJ ) ) Routput = true;
			//                                    if ( allOutput ) Routput = true;
												if ( Rdebug ) Routput = true;
												if ( intersectionOutput && Routput ) debugStream.Resume(); // enable logging
												debugStream << "      = R neighbour ====================" << std::endl;
												debugStream << "      basefunctions " << i << " " << j << std::endl;
												debugStream << "      faceQuadratureNeighbour.nop() " << faceQuadratureNeighbour.nop() << std::endl;
			#endif
												// sum over all quadrature points
												for ( size_t quad = 0; quad < faceQuadratureNeighbour.nop(); ++quad ) {
													// get x codim<0> and codim<1> coordinates
													const ElementCoordinateType xInside = faceQuadratureElement.point( quad );
													const ElementCoordinateType xOutside = faceQuadratureNeighbour.point( quad );
													const LocalIntersectionCoordinateType xLocal = faceQuadratureNeighbour.localPoint( quad );
													// get the integration factor
													const double elementVolume = intersectionGeoemtry.integrationElement( xLocal );
													// get the quadrature weight
													const double integrationWeight = faceQuadratureNeighbour.weight( quad );
													// compute \hat{u}_{p}^{P^{-}}(q_{j})\cdot n_{T}q_{i}
													const VelocityRangeType outerNormal = intersection.unitOuterNormal( xLocal );
													PressureRangeType q_j( 0.0 );
													pressureBaseFunctionSetNeighbour.evaluate( j, xOutside, q_j );
													PressureRangeType q_i( 0.0 );
													pressureBaseFunctionSetElement.evaluate( i, xInside, q_i );
													const double q_i_times_q_j = q_i * q_j;
													R_i_j += -1.0
														* D_11
														* elementVolume
														* integrationWeight
														* q_i_times_q_j;
			#ifndef NLOG
													debugStream << "      - quadPoint " << quad;
													Stuff::printFieldVector( xInside, "xInside", debugStream, "        " );
													Stuff::printFieldVector( xOutside, "xOutside", debugStream, "        " );
													Stuff::printFieldVector( xLocal, "xLocal", debugStream, "          " );
													Stuff::printFieldVector( outerNormal, "outerNormal", debugStream, "          " );
													debugStream << "\n        - elementVolume: " << elementVolume << std::endl;
													debugStream << "        - integrationWeight: " << integrationWeight;
													Stuff::printFieldVector( q_i, "q_i", debugStream, "        " );
													Stuff::printFieldVector( q_j, "q_j", debugStream, "        " );
													debugStream << "\n          - q_i_times_q_j: " << q_i_times_q_j << std::endl;
													debugStream << "          - R_" << i << "_" << j << "+=: " << R_i_j << std::endl;
			#endif
												} // done sum over all quadrature points
												// if small, should be zero
												if ( fabs( R_i_j ) < eps ) {
													R_i_j = 0.0;
												}
												else
													// add to matrix
													localRmatrixNeighbour.add( i, j, R_i_j );
			#ifndef NLOG
												Routput = false;
												debugStream.Suspend(); // disable logging
			#endif
											} // done computing R's neighbour surface integral
										} // done computing R's surface integrals
			//                        }
								} // done with those inside the grid

								// if we are on the boundary of the grid
								if ( !intersection.neighbor() && intersection.boundary() ) {


									//                                                                                                    // we will call this one
									// (H1)_{j} = \int_{\varepsilon\in\Epsilon_{D}^{T}}\hat{u}_{\sigma}^{RHS}()\cdot\tau_{j}\cdot n_{T}ds // H1's boundary integral
									if ( discreteModel_.hasVelocitySigmaFlux() ) {
										for ( int j = 0; j < numSigmaBaseFunctionsElement; ++j ) {
											double H1_j = 0.0;
			#ifndef NLOG
			//                                if ( j == logBaseJ ) H1output = true;
			//                                if ( allOutput ) H1output = true;
											if ( H1debug ) H1output = true;
											if ( intersectionOutput && H1output ) debugStream.Resume(); // enable logging
											debugStream << "      = H1 boundary ====================" << std::endl;
											debugStream << "      basefunction " << j << std::endl;
											debugStream << "      faceQuadratureElement.nop() " << faceQuadratureElement.nop() << std::endl;
			#endif
											// sum over all quadrature points
											for ( size_t quad = 0; quad < faceQuadratureElement.nop(); ++quad ) {
												// get x codim<0> and codim<1> coordinates
												const ElementCoordinateType x = faceQuadratureElement.point( quad );
												const VelocityRangeType xWorld = geometry.global( x );
												const LocalIntersectionCoordinateType xLocal = faceQuadratureElement.localPoint( quad );
												// get the integration factor
												const double elementVolume = intersectionGeoemtry.integrationElement( xLocal );
												// get the quadrature weight
												const double integrationWeight = faceQuadratureElement.weight( quad );
												// compute \hat{u}_{\sigma}^{RHS}()\cdot\tau_{j}\cdot n_{T}
												const VelocityRangeType outerNormal = intersection.unitOuterNormal( xLocal );
												SigmaRangeType tau_j( 0.0 );
												sigmaBaseFunctionSetElement.evaluate( j, x, tau_j );
												VelocityRangeType tau_j_times_normal( 0.0 );
												tau_j.mv( outerNormal, tau_j_times_normal );
												VelocityRangeType gD( 0.0 );
												discreteModel_.dirichletData( intersection, 0.0, xWorld,  gD );
												const double gD_times_tau_j_times_normal = gD * tau_j_times_normal;
												H1_j += elementVolume
													* integrationWeight
													* gD_times_tau_j_times_normal;
			#ifndef NLOG
												debugStream << "      - quadPoint " << quad;
												Stuff::printFieldVector( x, "x", debugStream, "        " );
												Stuff::printFieldVector( xLocal, "xLocal", debugStream, "        " );
												Stuff::printFieldVector( xWorld, "xWorld", debugStream, "        " );
												Stuff::printFieldVector( outerNormal, "outerNormal", debugStream, "        " );
												debugStream << "\n        - elementVolume: " << elementVolume << std::endl;
												debugStream << "        - integrationWeight: " << integrationWeight;
												Stuff::printFieldMatrix( tau_j, "tau_j", debugStream, "        " );
												Stuff::printFieldVector( gD, "gD", debugStream, "        " );
												Stuff::printFieldVector( tau_j_times_normal, "tau_j_times_normal", debugStream, "        " );
												debugStream << "\n        - gD_times_tau_j_times_normal: " << gD_times_tau_j_times_normal << std::endl;
												debugStream << "        - H1_" << j << "+=: " << H1_j << std::endl;
			#endif
											} // done sum over all quadrature points
											// if small, should be zero
											if ( fabs( H1_j ) < eps ) {
												H1_j = 0.0;
											}
											// add to rhs
											localH1rhs[ j ] += H1_j;
			#ifndef NLOG
											H1output = false;
											debugStream.Suspend(); // disable logging
			#endif
										} // done computing H1's boundary integral
									}

									//                                                                                                                   // we will call this one
									// (X)_{i,j} += \int_{\varepsilon\in\Epsilon_{D}^{T}}-\mu v_{i}\cdot\hat{\sigma}^{\sigma^{+}}(\tau_{j})\cdot n_{t}ds // X's boundary integral
									//                                                                                                                   // see also "X's volume integral", "X's element surface integral" and "X's neighbour surface integral" above
			#ifndef NLOG
									if ( Xdebug ) {
										debugStream.Resume(); // enable logging
										debugStream << "      = X boundary =====================" << std::endl;
										debugStream.Suspend();
									}
			#endif
			//                        if ( discreteModel_.hasSigmaFlux() ) {
										for ( int i = 0; i < numVelocityBaseFunctionsElement; ++i ) {
											for ( int j = 0; j < numSigmaBaseFunctionsElement; ++j ) {
												double X_i_j = 0.0;
			#ifndef NLOG
			//                                    if ( ( i == logBaseI ) && ( j == logBaseJ ) ) Xoutput = true;
			//                                    if ( allOutput ) Xoutput = true;
												if ( Xdebug ) Xoutput = true;
												if ( intersectionOutput && Xoutput ) debugStream.Resume(); // enable logging
												debugStream << "      = X boundary =====================" << std::endl;
												debugStream << "      basefunctions " << i << " " << j << std::endl;
												debugStream << "      faceQuadratureElement.nop() " << faceQuadratureElement.nop() << std::endl;
			#endif
												// sum over all quadrature points
												for ( size_t quad = 0; quad < faceQuadratureElement.nop(); ++quad ) {
													// get x codim<0> and codim<1> coordinates
													const ElementCoordinateType x = faceQuadratureElement.point( quad );
													const LocalIntersectionCoordinateType xLocal = faceQuadratureElement.localPoint( quad );
													// get the integration factor
													const double elementVolume = intersectionGeoemtry.integrationElement( xLocal );
													// get the quadrature weight
													const double integrationWeight = faceQuadratureElement.weight( quad );
													// compute -\mu v_{i}\cdot\hat{\sigma}^{\sigma^{+}}(\tau_{j})\cdot n_{t}
													const VelocityRangeType outerNormal = intersection.unitOuterNormal( xLocal );
													SigmaRangeType tau_j( 0.0 );
													sigmaBaseFunctionSetElement.evaluate( j, x, tau_j );
													VelocityRangeType v_i( 0.0 );
													velocityBaseFunctionSetElement.evaluate( i, x, v_i );
													VelocityRangeType tau_times_normal( 0.0 );
													tau_j.mv( outerNormal, tau_times_normal );
													const double v_i_times_tau_times_normal = v_i * tau_times_normal;
													X_i_j += -1.0
														* elementVolume
														* integrationWeight
														* mu
														* v_i_times_tau_times_normal;
			#ifndef NLOG
													debugStream << "      - quadPoint " << quad;
													Stuff::printFieldVector( x, "x", debugStream, "        " );
													Stuff::printFieldVector( xLocal, "xLocal", debugStream, "          " );
													Stuff::printFieldVector( outerNormal, "outerNormal", debugStream, "          " );
													debugStream << "\n        - elementVolume: " << elementVolume << std::endl;
													debugStream << "        - integrationWeight: " << integrationWeight;
													Stuff::printFieldVector( v_i, "v_i", debugStream, "        " );
													Stuff::printFieldMatrix( tau_j, "tau_j", debugStream, "        " );
													Stuff::printFieldVector( tau_times_normal, "tau_times_normal", debugStream, "        " );
													debugStream << "\n          - v_i_times_tau_times_normal: " << v_i_times_tau_times_normal << std::endl;
													debugStream << "          - X_" << i << "_" << j << "+=: " << X_i_j << std::endl;
			#endif
												} // done sum over all quadrature points
												// if small, should be zero
												if ( fabs( X_i_j ) < eps ) {
													X_i_j = 0.0;
												}
												else
													// add to matrix
													localXmatrixElement.add( i, j, X_i_j );
			#ifndef NLOG
												Xoutput = false;
												debugStream.Suspend(); // disable logging
			#endif
											}
										} // done computing X's boundary integral
			//                        }

									//                                                                                                           // we will call this one
									// (Y)_{i,j} += \int_{\varepsilon\in\Epsilon_{D}^{T}}-\mu v_{i}\cdot\hat{\sigma}^{U^{+}}(v_{j})\cdot n_{t}ds // Y's boundary integral
									//                                                                                                           // see also "Y's element surface integral" and "Y's neighbour surface integral" above
			#ifndef NLOG
									if ( Ydebug ) {
										debugStream.Resume(); // enable logging
										debugStream << "      = Y boundary =====================" << std::endl;
										debugStream.Suspend();
									}
			#endif
			//                        if ( discreteModel_.hasSigmaFlux() ) {
										for ( int i = 0; i < numVelocityBaseFunctionsElement; ++i ) {
											for ( int j = 0; j < numVelocityBaseFunctionsElement; ++j ) {
												double Y_i_j = 0.0;
			#ifndef NLOG
			//                                    if ( ( i == logBaseI ) && ( j == logBaseJ ) ) Youtput = true;
			//                                    if ( allOutput ) Youtput = true;
												if ( Ydebug ) Youtput = true;
												if ( intersectionOutput && Youtput ) debugStream.Resume(); // enable logging
												debugStream << "      = Y boundary =====================" << std::endl;
												debugStream << "      basefunctions " << i << " " << j << std::endl;
												debugStream << "      faceQuadratureElement.nop() " << faceQuadratureElement.nop() << std::endl;
			#endif
												// sum over all quadrature points
												for ( size_t quad = 0; quad < faceQuadratureElement.nop(); ++quad ) {
													// get x codim<0> and codim<1> coordinates
													const ElementCoordinateType x = faceQuadratureElement.point( quad );
													const ElementCoordinateType xInside = faceQuadratureElement.point( quad );
													const VelocityRangeType xWorld = geometry.global( x );
//													const ElementCoordinateType xOutside = faceQuadratureNeighbour.point( quad );
													const LocalIntersectionCoordinateType xLocal = faceQuadratureElement.localPoint( quad );
													// get the integration factor
													const double elementVolume = intersectionGeoemtry.integrationElement( xLocal );
													// get the quadrature weight
													const double integrationWeight = faceQuadratureElement.weight( quad );
													// compute -\mu v_{i}\cdot\hat{\sigma}^{U^{+}}(v_{j})\cdot n_{t}
													const VelocityRangeType outerNormal = intersection.unitOuterNormal( xLocal );
													VelocityRangeType v_j( 0.0 );
													velocityBaseFunctionSetElement.evaluate( j, x, v_j );
													VelocityRangeType v_i( 0.0 );
													velocityBaseFunctionSetElement.evaluate( i, x, v_i );
													const double v_i_times_v_j = v_i * v_j;
													Y_i_j += C_11
														* elementVolume
														* integrationWeight
														* mu
														* v_i_times_v_j;
													//calc (\beta * n) *  ( \hat{u}^c_h * v )
													//boundary edge
													VelocityRangeType beta_eval;
													beta_.evaluate( xWorld, beta_eval );
													const double beta_times_normal = beta_eval * outerNormal;
													VelocityRangeType flux_value;
													if ( beta_times_normal < 0 ) {
//														velocityBaseFunctionSetElement.evaluate( i, xInside, flux_value );
														VelocityRangeType gD( 0.0 );
														discreteModel_.dirichletData( intersection, 0.0, xWorld,  gD );
														flux_value = gD;
													}
													else {
														velocityBaseFunctionSetElement.evaluate( i, xInside, flux_value );
//														VelocityRangeType gD( 0.0 );
//														discreteModel_.dirichletData( intersection, 0.0, xWorld,  gD );
//														flux_value = gD;
													}
													const double flux_times_v_i = flux_value * v_i;
													Y_i_j += elementVolume
														* integrationWeight
														* beta_times_normal
														* flux_times_v_i;
			#ifndef NLOG
													debugStream << "      - quadPoint " << quad;
													Stuff::printFieldVector( x, "x", debugStream, "        " );
													Stuff::printFieldVector( xLocal, "xLocal", debugStream, "          " );
													Stuff::printFieldVector( outerNormal, "outerNormal", debugStream, "          " );
													debugStream << "\n        - elementVolume: " << elementVolume << std::endl;
													debugStream << "        - integrationWeight: " << integrationWeight;
													Stuff::printFieldVector( v_i, "v_i", debugStream, "        " );
													Stuff::printFieldVector( v_j, "v_j", debugStream, "        " );
													debugStream << "\n          - v_i_times_v_j: " << v_i_times_v_j << std::endl;
													debugStream << "          - Y_" << i << "_" << j << "+=: " << Y_i_j << std::endl;
			#endif
												} // done sum over all quadrature points
												// if small, should be zero
												if ( fabs( Y_i_j ) < eps ) {
													Y_i_j = 0.0;
												}
												else
													// add to matrix
													localYmatrixElement.add( i, j, Y_i_j );
			#ifndef NLOG
												Youtput = false;
												debugStream.Suspend(); // disable logging
			#endif
											}
										} // done computing Y's boundary integral
			//                        }

									//                                                                                                  // we will call this one
									// (Z)_{i,j} += \int_{\varepsilon\in\Epsilon_{D}^{T}}\hat{p}^{P^{+}}(q_{j})\cdot v_{i}\cdot n_{T}ds // Z's boundary integral
									//                                                                                                  // see also "Z's volume integral", "Z's element surface integral" and "Z's neighbour surface integral" above
			#ifndef NLOG
									if ( Zdebug ) {
										debugStream.Resume(); // enable logging
										debugStream << "      = Z boundary =====================" << std::endl;
										debugStream.Suspend();
									}
			#endif
			//                        if ( discreteModel_.hasPressureFlux() ) {
										for ( int i = 0; i < numVelocityBaseFunctionsElement; ++i ) {
											// compute the boundary integral
											for ( int j = 0; j < numPressureBaseFunctionsElement; ++j ) {
												double Z_i_j = 0.0;
			#ifndef NLOG
			//                                    if ( ( i == logBaseI ) && ( j == logBaseJ ) ) Zoutput = true;
			//                                    if ( allOutput ) Zoutput = true;
												if ( Zdebug ) Zoutput = true;
												if ( intersectionOutput && Zoutput ) debugStream.Resume(); // enable logging
												debugStream << "      = Z boundary =====================" << std::endl;
												debugStream << "      basefunctions " << i << " " << j << std::endl;
												debugStream << "      faceQuadratureElement.nop() " << faceQuadratureElement.nop() << std::endl;
			#endif
												// sum over all quadrature points
												for ( size_t quad = 0; quad < faceQuadratureElement.nop(); ++quad ) {
													// get x codim<0> and codim<1> coordinates
													const ElementCoordinateType x = faceQuadratureElement.point( quad );
													const LocalIntersectionCoordinateType xLocal = faceQuadratureElement.localPoint( quad );
													// get the integration factor
													const double elementVolume = intersectionGeoemtry.integrationElement( xLocal );
													// get the quadrature weight
													const double integrationWeight = faceQuadratureElement.weight( quad );
													// compute \hat{p}^{P^{+}}(q_{j})\cdot v_{i}\cdot n_{T}
													const VelocityRangeType outerNormal = intersection.unitOuterNormal( xLocal );
													VelocityRangeType v_i( 0.0 );
													velocityBaseFunctionSetElement.evaluate( i, x, v_i );
													PressureRangeType q_j( 0.0 );
													pressureBaseFunctionSetElement.evaluate( j, x, q_j );
													const double v_i_times_normal = v_i * outerNormal;
													const double q_j_times_v_times_normal = q_j * v_i_times_normal;
													Z_i_j += elementVolume
														* integrationWeight
														* q_j_times_v_times_normal;
			#ifndef NLOG
													debugStream << "      - quadPoint " << quad;
													Stuff::printFieldVector( x, "x", debugStream, "        " );
													Stuff::printFieldVector( xLocal, "xLocal", debugStream, "          " );
													Stuff::printFieldVector( outerNormal, "outerNormal", debugStream, "          " );
													debugStream << "\n        - elementVolume: " << elementVolume << std::endl;
													debugStream << "        - integrationWeight: " << integrationWeight;
													Stuff::printFieldVector( v_i, "v_i", debugStream, "        " );
													Stuff::printFieldVector( q_j, "q_j", debugStream, "        " );
													debugStream << "\n          - q_j_times_v_times_normal: " << q_j_times_v_times_normal << std::endl;
													debugStream << "          - Z_" << i << "_" << j << "+=: " << Z_i_j << std::endl;
			#endif
												} // done sum over all quadrature points
												// if small, should be zero
												if ( fabs( Z_i_j ) < eps ) {
													Z_i_j = 0.0;
												}
												else
													// add to matrix
													localZmatrixElement.add( i, j, Z_i_j );
			#ifndef NLOG
												Zoutput = false;
												debugStream.Suspend(); // disable logging
			#endif
											}
										} // done computing Z's boundary integral
			//                        }

									//                                                                                                                 // we will call this one
									// (H2)_{j} += \int_{\varepsilon\in\Epsilon_{D}^{T}}\left( \mu v_{j}\cdot\hat{\sigma}^{RHS}()\cdot n_{T}ds         // H2's 1st boundary integral
									//                                                         -\hat{p}^{RHS}()\cdot v_{j}\cdot n_{T}ds        \right) // H2's 2nd boundary integral
									//                                                                                                                 // see also "H2's volume integral" above
									if ( discreteModel_.hasSigmaFlux() && discreteModel_.hasPressureFlux() ) {
										for ( int j = 0; j < numVelocityBaseFunctionsElement; ++j ) {
											double H2_j = 0.0;
			#ifndef NLOG
			//                                if ( j == logBaseJ ) H2output = true;
			//                                if ( allOutput ) H2output = true;
											if ( H2debug ) H2output = true;
											if ( intersectionOutput && H2output ) debugStream.Resume(); // enable logging
											debugStream << "      = H2 boundary ====================" << std::endl;
											debugStream << "      basefunction " << j << std::endl;
											debugStream << "      faceQuadratureElement.nop() " << faceQuadratureElement.nop() << std::endl;
			#endif
											// sum over all quadrature points
											for ( size_t quad = 0; quad < faceQuadratureElement.nop(); ++quad ) {
												// get x codim<0> and codim<1> coordinates
												const ElementCoordinateType x = faceQuadratureElement.point( quad );
												const LocalIntersectionCoordinateType xLocal = faceQuadratureElement.localPoint( quad );
												const VelocityRangeType globalX = geometry.global( x );
												// get the integration factor
												const double elementVolume = intersectionGeoemtry.integrationElement( xLocal );
												// get the quadrature weight
												const double integrationWeight = faceQuadratureElement.weight( quad );
												// prepare
												const VelocityRangeType outerNormal = intersection.unitOuterNormal( xLocal );
												VelocityRangeType v_j( 0.0 );
												velocityBaseFunctionSetElement.evaluate( j, x, v_j );
												// compute \mu v_{j}\cdot\hat{\sigma}^{RHS}()\cdot n_{T}
			//                                    if ( discreteModel_.hasSigmaFlux() ) {
													const VelocityRangeType xIntersectionGlobal = intersection.intersectionSelfLocal().global( xLocal );
													const VelocityRangeType xWorld = geometry.global( xIntersectionGlobal );
													VelocityRangeType gD( 0.0 );
													discreteModel_.dirichletData( intersection, 0.0, xWorld, gD );
													SigmaRangeType gD_times_normal( 0.0 );
													gD_times_normal = dyadicProduct( gD, outerNormal );
													VelocityRangeType gD_times_normal_times_normal( 0.0 );
													gD_times_normal.mv( outerNormal, gD_times_normal_times_normal );
													const double v_j_times_gD_times_normal_times_normal= v_j * gD_times_normal_times_normal;
													H2_j += C_11
														* elementVolume
														* integrationWeight
														* mu
														* v_j_times_gD_times_normal_times_normal;
			#ifndef NLOG
													debugStream << "      - quadPoint " << quad;
													Stuff::printFieldVector( x, "x", debugStream, "        " );
													Stuff::printFieldVector( xLocal, "xLocal", debugStream, "        " );
													Stuff::printFieldVector( globalX, "globalX", debugStream, "        " );
													Stuff::printFieldVector( outerNormal, "outerNormal", debugStream, "        " );
													debugStream << "\n        - elementVolume: " << elementVolume << std::endl;
													debugStream << "        - integrationWeight: " << integrationWeight;
													Stuff::printFieldVector( v_j, "v_j", debugStream, "        " );
													Stuff::printFieldVector( gD, "gD", debugStream, "        " );
													Stuff::printFieldMatrix( gD_times_normal, "gD_times_normal", debugStream, "        " );
													Stuff::printFieldVector( gD_times_normal_times_normal, "gD_times_normal_times_normal", debugStream, "        " );
													debugStream << "\n        - v_j_times_gD_times_normal_times_normal: " << v_j_times_gD_times_normal_times_normal << std::endl;
													debugStream << "        - H2_" << j << "+=: " << H2_j;
			#endif
			//                                    }
												// done computing H2's 1st boundary integral
												// compute -\hat{p}^{RHS}()\cdot v_{j}\cdot n_{T}
			//                                    if ( discreteModel_.hasPressureFlux() ) {
													const double v_j_times_normal = v_j * outerNormal;
													const double flux_times_v_j_times_n_t = 0.0 * v_j_times_normal;
													H2_j += -1.0
														* elementVolume
														* integrationWeight
														* flux_times_v_j_times_n_t;
			#ifndef NLOG
													debugStream << "        - v_j_times_normal: " << v_j_times_normal << std::endl;
													debugStream << "        - flux_times_v_j_times_n_t: " << flux_times_v_j_times_n_t << std::endl;
													debugStream << "        - H2_" << j << "+=: " << H2_j << std::endl;
			#endif
			//                                    }
												// done computing H2's 2nd boundary integral
											} // done sum over all quadrature points
											// if small, should be zero
											if ( fabs( H2_j ) < eps ) {
												H2_j = 0.0;
											}
											// add to rhs
											localH2rhs[ j ] += H2_j;
			#ifndef NLOG
											H2output = false;
											debugStream.Suspend(); // disable logging
			#endif
										} // done computing H2's boundary integrals
									}

									//                                                                                        // we will call this one
									// (H3)_{j} = \int_{\varepsilon\in\Epsilon_{D}^{T}}-\hat{u}_{p}^{RHS}()\cdot n_{T}q_{j}ds // H3's boundary integral
			//                        if ( discreteModel_.hasVelocityPressureFlux() ) {
										for ( int j = 0; j < numPressureBaseFunctionsElement; ++j ) {
											double H3_j = 0.0;
			#ifndef NLOG
			//                                if ( j == logBaseJ ) H3output = true;
			//                                if ( allOutput ) H3output = true;
											if ( H3debug ) H3output = true;
											if ( intersectionOutput && H3output ) debugStream.Resume(); // enable logging
											debugStream << "      = H3 boundary ====================" << std::endl;
											debugStream << "      basefunction " << j << std::endl;
											debugStream << "      faceQuadratureElement.nop() " << faceQuadratureElement.nop() << std::endl;
			#endif
											// sum over all quadrature points
											for ( size_t quad = 0; quad < faceQuadratureElement.nop(); ++quad ) {
												// get x codim<0> and codim<1> coordinates
												const ElementCoordinateType x = faceQuadratureElement.point( quad );
												const LocalIntersectionCoordinateType xLocal = faceQuadratureElement.localPoint( quad );
												const VelocityRangeType xWorld = geometry.global( x );
												// get the integration factor
												const double elementVolume = intersectionGeoemtry.integrationElement( xLocal );
												// get the quadrature weight
												const double integrationWeight = faceQuadratureElement.weight( quad );
												// compute -\hat{u}_{p}^{RHS}()\cdot n_{T}q_{j}
												const VelocityRangeType outerNormal = intersection.unitOuterNormal( xLocal );
												VelocityRangeType gD( 0.0 );
												discreteModel_.dirichletData( intersection, 0.0, xWorld, gD );
												const double gD_times_normal = gD * outerNormal;
												PressureRangeType q_j( 0.0 );
												pressureBaseFunctionSetElement.evaluate( j, x, q_j );
												const double q_j_times_gD_times_normal = q_j * gD_times_normal;
												H3_j += elementVolume
													* integrationWeight
													* q_j_times_gD_times_normal;
			#ifndef NLOG
												debugStream << "      - quadPoint " << quad;
												Stuff::printFieldVector( x, "x", debugStream, "        " );
												Stuff::printFieldVector( xLocal, "xLocal", debugStream, "        " );
												Stuff::printFieldVector( outerNormal, "outerNormal", debugStream, "        " );
												debugStream << "\n        - elementVolume: " << elementVolume << std::endl;
												debugStream << "        - integrationWeight: " << integrationWeight;
												Stuff::printFieldVector( q_j, "q_j", debugStream, "        " );
												Stuff::printFieldVector( gD, "gD", debugStream, "        " );
												debugStream << "\n        - gD_times_normal: " << gD_times_normal << std::endl;
												debugStream << "        - q_j_times_gD_times_normal: " << q_j_times_gD_times_normal << std::endl;
												debugStream << "        - H3_" << j << "+=: " << H3_j << std::endl;
			#endif
											} // done sum over all quadrature points
											// if small, should be zero
											if ( fabs( H3_j ) < eps ) {
												H3_j = 0.0;
											}
											// add to rhs
											localH3rhs[ j ] += H3_j;
			#ifndef NLOG
											H3output = false;
											debugStream.Suspend(); // disable logging
			#endif
										} // done computing H3's boundary integral
			//                        }

								} // done with those on the boundary
			#ifndef NLOG
								if ( intersectionOutput ) debugStream.Resume(); // enable logging
								debugStream << "    - ==== done calculations on intersection " << intersectionNR << std::endl;
								debugStream.Suspend(); // disable logging
								intersectionOutput = false;
								++intersectionNR;
			#endif
							} // done walking the neighbours
			//#endif //no_surface_ints

			#ifndef NLOG
							intersectionNR = 0;
							if ( entityOutput ) debugStream.Resume(); // enable logging
							debugStream << "  - == done calculations on entity " << entityNR << std::endl;
							debugStream.Suspend(); // disable logging
							entityOutput = false;
							++entityNR;
			#endif
						} // done walking the grid


			#ifndef NLOG
						infoStream.Resume();
						if ( numberOfEntities > 19 ) {
							infoStream << "]";
						}
						infoStream << "\n- gridwalk done" << std::endl << std::endl;
						infoStream.Suspend();

			#endif
						// do profiling
						profiler().StopTiming("Pass -- ASSEMBLE");

						profiler().StartTiming("Pass -- SOLVER");

						// do solving
			#ifndef NLOG
						infoStream.Resume();
						infoStream << "Solving system with " << dest.discreteVelocity().size() << " + " << dest.discretePressure().size() << " unknowns" << std::endl;
			#endif
						InvOpType op;
			#ifdef USE_NESTED_CG_SOLVER
						AltInvOpType m_op;
						if ( Parameters().getParam( "use_nested_cg_solver", false ) )
							info_ = m_op.solve( arg, dest, Xmatrix, MInversMatrix, Ymatrix, Ematrix, Rmatrix, Zmatrix, Wmatrix, H1rhs, H2rhs, H3rhs );
						else
			#endif
						info_ = op.solve( arg, dest, Xmatrix, MInversMatrix, Ymatrix, Ematrix, Rmatrix, Zmatrix, Wmatrix, H1rhs, H2rhs, H3rhs );

						// do profiling
						profiler().StopTiming("Pass -- SOLVER");
					} // end of apply

					virtual void compute( const TotalArgumentType& /*arg*/, DestinationType& /*dest*/ ) const
					{}

					virtual void allocateLocalMemory()
					{}

			#ifdef HAS_RUN_INFO
					void getRuninfo( RunInfo& info )
					{
						info.iterations_inner_avg = int( info_.iterations_inner_avg );
						info.iterations_inner_min = info_.iterations_inner_min;
						info.iterations_inner_max = info_.iterations_inner_max;
						info.iterations_outer_total = info_.iterations_outer_total;
						info.max_inner_accuracy = info_.max_inner_accuracy;
					}
			#endif

				private:
					DiscreteModelType& discreteModel_;
					GridPartType& gridPart_;
			//        NonPeriodicGridPartType& nonPeriodicGridPart_;
					DiscreteStokesFunctionSpaceWrapperType& spaceWrapper_;
					DiscreteVelocityFunctionSpaceType& velocitySpace_;
					DiscretePressureFunctionSpaceType& pressureSpace_;
					DiscreteSigmaFunctionSpaceType sigmaSpace_;
					const ExtraDataFunctionType& beta_;
					mutable SaddlepointInverseOperatorInfo info_;

					/**
					 *  \todo   doc
					 **/

					/**
					 *  \brief  dyadic product
					 *
					 *          Implements \f$\left(arg_{1} \otimes arg_{2}\right)_{i,j}:={arg_{1}}_{i} {arg_{2}}_{j}\f$
					 **/
					static SigmaRangeType dyadicProduct(   const VelocityRangeType& arg1,
													const VelocityRangeType& arg2 )
					{
						SigmaRangeType ret( 0.0 );
						typedef typename SigmaRangeType::RowIterator
							MatrixRowIteratorType;
						typedef typename VelocityRangeType::ConstIterator
							ConstVectorIteratorType;
						typedef typename VelocityRangeType::Iterator
							VectorIteratorType;
						MatrixRowIteratorType rItEnd = ret.end();
						ConstVectorIteratorType arg1It = arg1.begin();
						for ( MatrixRowIteratorType rIt = ret.begin(); rIt != rItEnd; ++rIt ) {
							ConstVectorIteratorType arg2It = arg2.begin();
							VectorIteratorType vItEnd = rIt->end();
							for (   VectorIteratorType vIt = rIt->begin();
									vIt != vItEnd;
									++vIt ) {
								*vIt = *arg1It * *arg2It;
								++arg2It;
							}
							++arg1It;
						}
						return ret;
					}

//					template <class
//					static divergence_of_dyadic_product
//
					// VelocityRangeType is expected to be a FieldVector,
					// SigmaJacobianRangeType to be a Matrixmapping and
					// SigmaJacobianRangeType[i] to be a FieldVector
					//! \todo   doc me
					static SigmaJacobianRangeType prepareVelocityRangeTypeForSigmaDivergence( const VelocityRangeType& arg )
					{
						SigmaJacobianRangeType ret( 0.0 );
						assert( arg.dim() == ret[0].dim() );
						for ( int i = 0; i < int(arg.dim()) ; ++i ) {
							for ( int j = 0; j < int(arg.dim()); ++j ) {
								VelocityRangeType row( 0.0 );
								row[ j ] = arg[ i ];
								ret[ i * arg.dim() + j ] = row;
							}
						}
						return ret;
					}

					//! \todo   doc me
					static VelocityJacobianRangeType preparePressureRangeTypeForVelocityDivergence( const PressureRangeType& arg )
					{
						VelocityJacobianRangeType ret( 0.0 );
						for ( unsigned int i = 0; i < ret[0].dim(); ++i ) {
							VelocityRangeType row( 0.0 );
							row[ i ] = arg;
							ret[ i ] = row;
						}
						return ret;
					}


			};
		}//end namespace NonlinearStep
	}//end namespace NavierStokes
}//end namespace Dune

#endif  // end of stokespass.hh
