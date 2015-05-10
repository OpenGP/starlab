*****************************************************
****************** PQP++ Nov. 2010 ******************
*****************************************************

-- Introduction --
------------------
PQP++ is a port of PQP v1.3, which was originally developed by S. Gottschalk and E. Larsen at University of North Carolina at Chapel Hill.

-- Changes --
-------------
Following changes have been made:
* All classes and enums are within the namespace PQP
* All global methods are put in classes in order to allow multiple instanciations.

-- Usage --
-----------
Instead of calling the global functions ::PQP_Collide() and ::PQP_Distance() create an instance of the class PQP_Checker and use the according methods of this instance:
  PQP::PQP_Checker *pqpChecker = new PQP_Checker();
  PQP::PQP_CollideResult result;
	pqpChecker->PQP_Collide(&result, R1, T1, m1, R2, T2, m2, PQP::PQP_FIRST_CONTACT);
or
  PQP::PQP_DistanceResult pqpResult;
  pqpChecker->PQP_Distance(&pqpResult, Rotation1, Translation1, m1, Rotation2, Translation2, m2, 0, 0);
  
  
  
-- Contact Information --
-------------------------
Nikolaus Vahrenkamp
Institute for Anthropomatics
HIS - Chair Prof. Ruediger Dillmann
Karlsruhe Institute of Technology (KIT)
Germany

mail: vahrenkamp@kit.edu

