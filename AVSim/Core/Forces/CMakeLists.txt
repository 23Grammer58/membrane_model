set(HEPOT_COLLECTION_SOURCES "")
set(HEPOT_COLLECTION_HEADERS NeoGook.h Gent.h SVKirchhoff.h MooneyRivlin.h MayYin.h HolzapfelOgdenGasser.h HGO_GSTModel.h HyperElasticModel.h)
set(FORCE_SOURCES Forces.cpp ForcesCommon.cpp
                Integrator1D.cpp HyperElastic.cpp BendHyperElastic.cpp
                MSM.cpp LinearElastic.cpp ${HEPOT_COLLECTION_SOURCES}
                Pressure.cpp EdgeLoad.cpp SpurnPlane.cpp ClampedBndPenalty.cpp VertexLoad.cpp
                SDFForce.cpp
                HyperElasticResponse.cpp
                PARENT_SCOPE)
set(CONTACT_FORCE_SOURCES ContactForce.cpp PARENT_SCOPE)
set(CONTACT_FORCE_HEADERS ContactForce.h PARENT_SCOPE)                
set(FORCE_HEADERS Forces.h ForcesCommon.h
                Integrator1D.h HyperElastic.h HyperElastic.inl BendHyperElastic.h MSM.h LinearElastic.h
                ${HEPOT_COLLECTION_HEADERS}
                Pressure.h Pressure.inl EdgeLoad.h SpurnPlane.h ClampedBndPenalty.h VertexLoad.h VertDampForce.h
                SDFForce.h
                HyperElasticResponse.h MaterialThicknessHelper.h
                PARENT_SCOPE)
set(DATA_DRIVEN_SOURCES DataDriven/MaterialResponse.cpp 
                        DataDriven/ConjugatePairSearcher.cpp 
                        DataDriven/LinearElasticInterpolant.cpp
                        DataDriven.cpp #DataGatherer.cpp
                        PARENT_SCOPE) 
set(DATA_DRIVEN_HEADERS DataDriven.h #DataGatherer.h
                        DataDriven/MaterialResponse.h 
                        DataDriven/ConjugatePairSearcher.h 
                        DataDriven/LinearElasticInterpolant.h 
                        DataDriven/InterpolatorBase.h 
                        DataDriven/KNearestInterpolant.h DataDriven/KNearestInterpolant.inl
                        PARENT_SCOPE)                                       


