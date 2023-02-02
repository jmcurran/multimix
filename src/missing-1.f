	PROGRAM MISSING

*	This program fits a mixture of multivariate distributions
*	using the EM algorithm to a data set that has data missing at
*	random. The data file contains both categorical and continuous
*	variables. The missing data must be coded to the value set in
*	NOTE (2).

*	If the program does not converge after iter=200 iterations,
*	the estimates of the parameters will be entered into
*	EMPARAMEST.OUT. This file can then be used as the parameter
*	input file for PROGRAM MISSING if desired.

*	The assignment of the observations to groups (IGP(i)) and the
*	posterior probabilities (Zij's) are entered into GROUPS.OUT.
*	This file can be used in MINITAB etc. for further analysis.

*	NOTE:

*	(1) THIS PROGRAM REQUIRES VARIABLES IN A PARTITION TO BE STORED
*	CONTIGUOUSLY. HENCE THE DATA IS READ IN WITH THE VARIABLE ORDER
*	BEING SPECIFIED BY JP(J). IVARTYPE(J) AND NCAT(J) BOTH REFER TO
*	THE REARRANGED DATA

*	(2) The code for missing data is -999. imiss=-999.

*	(3) IF SDENS=0, THEN Z(II,K) IS SET TO 0.01

*	(4) THE PROGRAM CURRENTLY HAS A MAXIMUM OF
*		1500 OBSERVATIONS 		(iob=1500)
*		6 GROUPS	  		(ik6=6)
*		15 ATTRIBUTES & 15 PARTITIONS 	(ip15=15)
*		10 LEVELS OF CATEGORIES 	(im10=10)
*		200 ITERATIONS FOR CONVERGENCE 	(iter=200)
*         ******ALTER IF REQUIRED******
*	  (REMEMBER TO ALTER PARAMETERS IN DETINV ALSO)

*	The parameter file contains:-
*	NG - the number of groups
*	NOBS - the number of observations
*	NVAR - the number of variables
*	NPAR - the number of partitions
*	ISPEC - an indicator variable for a specified grouping of the
*		observations (1 = observations are not specified into
*		groups, 2 = observations are specified into groups)
*	JP(j) - column in the data array in which the jTH variable of
*	         the file will be stored
*	IP(l) - number of variables in the lTH partition, l=1,NPAR
*	IPC(l) - number of continuous variables in partition l, l=1,NPAR
*	ISV(l) - indicator starting value for the partition, l=1,NPAR
*	IEV(l) - indicator end value for the partition, l=1,NPAR
*	IPARTYPE(l) - indicator giving the type of model each partition is
*		(1 = categorical, 2 = MVN, 3 = location models), l=1,NPAR
*	IVARTYPE(j) - an indicator variable for type of variable, j=1,NVAR
*		(1 = categorical variable, 2 = continuous variable,
*		3 = categorical variable involved in location model,
*		4 = continuous variable involved in the location model)
*	NCAT(j) - number of categories of jTH variable.
*		(For continuous variables, set NCAT(J)=0)
*	PI(k) - estimated mixing proportions for each group, k=1,NG
*	THETA(K,J,M) - estimated probability that the jTH categorical
*			variable is at level M, given that in group k
*	EMUL(k,l,j,m)- estimated mean vector for each group, each
*	partition, and each location of the location model variables
*	EMU(K,L,J) -estimated mean vector for each group and each
*	partition of continuous variables
*	VARIX(K,L,I,J)   -estimated covariance matrices for each group

 	IMPLICIT DOUBLE PRECISION(A-H,O-Z)
	CHARACTER*16 infile, datafile
	PARAMETER (PIE=3.141592653589792, iob=1500, ip15=15, ik6=6,
     :           im10=10, iter=200, imiss=-9)
	DIMENSION EMU(ik6,ip15,ip15),VAR(ik6,ip15,ip15,ip15),
     :  Z(iob,ik6), PI(ik6), DENS(ik6,ip15), VARIN(ik6,ip15,ip15,ip15),
     :  ZSUM(ik6), XSUM(ik6,ip15), ADET(ik6,ip15), CLOGLI(iter),
     :  VARIX(ik6,ip15,ip15,ip15), APRODENS(ik6), IP(ip15), IPC(ip15),
     :  ISV(ip15), IEV(ip15), IPARTYPE(ip15), EMUL(ik6,ip15,ip15,im10),
     :  THETA(ik6,ip15,0:im10), PROB(ik6,iob),IVARTYPE(ip15),NCAT(ip15),
     :  IM(ip15), IGP(iob), NUM(ik6), XSUM2(ik6,ip15,ip15),ZM(ik6,ip15),
     :  IGRP(iob), JP(ip15), IX(iob,ip15), X(iob,ip15),
     :  IMISSCODE(iob,ip15), IIX(iob,ip15,im10), ISWPCOL(ip15),
     :	AUGCOV(ik6,ip15,0:ip15,0:ip15),EN(ik6,ip15,0:ip15,0:ip15),
     :  VAROBS(ik6,ip15,ip15,ip15), ESTX(ik6,iob,ip15),
     :  THETA2(ik6,ip15,im10)
	PRINT *, ' MIXTURE ESTIMATION BY EM'
	PRINT *, '-------------------------------'
	PRINT *, 'Data file: '
	READ (*,10) datafile
10	FORMAT (A16)
	PRINT *, 'Parameter file: '
	READ (*,10) infile
	OPEN (7, FILE='GENERAL.OUT', STATUS='NEW')
	OPEN (8, FILE=datafile, STATUS='OLD')
	OPEN (9, FILE=infile, STATUS='OLD')
	OPEN(12, FILE='GROUPS.OUT',STATUS='NEW')
	READ (9,*) NG, NOBS, NVAR, NPAR,ISPEC
	READ (9,*) (JP(J),J=1,NVAR)
	READ (9,*) (IP(L),L=1,NPAR)
	READ (9,*) (IPC(L), L=1,NPAR)
	READ (9,*) (ISV(L),L=1,NPAR)
	READ (9,*) (IEV(L),L=1,NPAR)
	READ (9,*) (IPARTYPE(L),L=1,NPAR)
	READ (9,*) (IVARTYPE(J),J=1,NVAR)
	READ (9,*) (NCAT(J),J=1,NVAR)

*	read in estimates of the parameters if grouping not specified

	IF(ISPEC.EQ.1) THEN
	  READ(9,*) (PI(K),K=1,NG)
	  DO K=1,NG	!12
	  DO L=1,NPAR	!12
		IF (IPARTYPE(L).EQ.1) THEN
			DO J=ISV(L),IEV(L)	!13
			  READ(9,*) (THETA(K,J,M),M=1,NCAT(J))
  			END DO	!13

		ELSE IF(IPARTYPE(L).EQ.2) THEN
			READ(9,*)(EMU(K,L,J),J=1,IPC(L))
		ELSE IF(IPARTYPE(L).EQ.3) THEN
			DO J=ISV(L),IEV(L)	!14
			  IF(IVARTYPE(J).EQ.3) THEN
				READ(9,*)(THETA(K,J,M),M=1,NCAT(J))
				IM(L)=NCAT(J)
			   END IF
  			END DO	!14

			J1=0
			DO J=ISV(L),IEV(L)	!15
			  IF(IVARTYPE(J).EQ.4) THEN
			   J1=J1+1
			   READ(9,*)
     :			    (EMUL(K,L,J1,M),M=1,IM(L))
			  END IF
  			END DO	!15

		END IF
  	END DO	!12
  	END DO	!12

	DO K=1,NG	!16
	DO L=1,NPAR	!16
		IF(IPARTYPE(L).NE.1)
     :		READ(9,*)((VARIX(K,L,I,J),J=1,IPC(L)),I=1,IPC(L))
    	END DO	!16
    	END DO	!16

	ELSE
	  READ(9,*) (IGRP(I),I=1,NOBS)
	  DO I=1,NOBS	!18
	  DO K=1,NG	!18
	      Z(I,K)=0.0
  	  END DO	!18
  	  END DO	!18
	  DO I=1,NOBS	!31
	    IK=IGRP(I)
	    Z(I,IK)=1.0
  	  END DO	!31

	  DO L=1,NPAR	!33
	     IF(IPARTYPE(L).EQ.3) THEN
		DO J=ISV(L),IEV(L)	!32
		  IF (IVARTYPE(J).EQ.3) IM(L)=NCAT(J)
  		END DO	!32

	    END IF
  	  END DO	!33

	END IF

*	read in the observations

	DO I=1,NOBS	!17
		READ(8,* ) (X(I,JP(J)),J=1,NVAR)
  	END DO	!17


*	Check for the presence of missing values, and create
*	IMISSCODE(I,J), where IMISSCODE(I,J)=1 if observation i variable j is
*	missing, and IMISSCODE(I,J)=0 otherwise

	DO I=1,NOBS	!34
	DO J=1,NVAR	!34
		IVAL=NINT(X(I,J))
		IF (IVAL.EQ.IMISS) THEN
			IMISSCODE(I,J)=1
		ELSE
			IMISSCODE(I,J)=0
		END IF
  	END DO	!34
  	END DO	!34

*	Check the categorical variable to see if 0 is a category.
*	Add 1 to the variables if it is so that Xij=0 is taken as
*	level 1, Xij=1 is level 2 etc.

	DO J=1,NVAR	!80
	    IF((IVARTYPE(J).EQ.1).OR.(IVARTYPE(J).EQ.3)) THEN
		IMIN=NINT(X(1,J))
		IF (IMIN.NE.0) THEN
		   DO I=2,NOBS	!81
			INEXT=NINT(X(I,J))
			IF(INEXT.EQ.0) IMIN=INEXT
  		   END DO	!81

		END IF
		IF (IMIN.EQ.0) THEN
		    DO I=1,NOBS	!82
			IF (IMISSCODE(I,J).EQ.0) THEN
		 	   IX(I,J)=NINT(X(I,J)) + 1
			ELSE
			   IX(I,J)=NINT(X(I,J))
			END IF
  		    END DO	!82

		ELSE
		    DO I=1,NOBS	!83
			IX(I,J)=NINT(X(I,J))
  		    END DO	!83

		END IF
	   END IF
  	END DO	!80


*	For the categorical variables, create the matrix IIX(I,J,M) where
*	x(i,j,m)=1 if obs. i, variable j is at level m,
*	x(i,j,m)=0 if not at level m.

	DO I=1,NOBS	!84
	DO J=1,NVAR	!84
	    IF((IVARTYPE(J).EQ.1).OR.(IVARTYPE(J).EQ.3)) THEN
		DO M=1,NCAT(J)	!85
		     IF (IX(I,J).EQ.M) THEN
			IIX(I,J,M)=1
		     ELSE IF (IX(I,J).EQ.IMISS) THEN
			IIX(I,J,M)=9
		     ELSE
			IIX(I,J,M)=0
		     END IF
  		END DO	!85

	     END IF
  	END DO	!84
  	END DO	!84

*	print the input values

	WRITE(7,50)NG,NOBS,NVAR,NPAR
50	FORMAT(/,' NO OF GROUPS IS ',I3,/,' NO OF OBSERVATIONS IS ',
     :  I5,/,' NO OF VARIABLES ',I3,/,' NO OF PARTITIONS', I3)
	WRITE(7,51)(IP(L),L=1,NPAR)
51	FORMAT(/,' THE NO OF VARIABLES IN EACH PARTITION IS',/,10I3)

*	Send to the M step if ISPEC = 2 (groups specified)
	IF(ISPEC.EQ.2) GO TO 100
*	otherwise, print the parameter estimates
	WRITE(7,52)
52	FORMAT(/,' MIXING PROPORTIONS')
	WRITE(7,53)(PI(K),K=1,NG)
53	FORMAT(/,10F6.3)
	DO K=1,NG	!56
	DO L=1,NPAR	!56
		IF(IPARTYPE(L).EQ.1) THEN
		  WRITE(7,54) K
54		  FORMAT(/,2X,' THETA(K,J,M) FOR GROUP ',I2)
		  DO J=ISV(L),IEV(L)	!64
		      WRITE(7,55)(THETA(K,J,M),M=1,NCAT(J))
55		      FORMAT(10F8.4)
  		  END DO	!64

		ELSE IF(IPARTYPE(L).EQ.2) THEN
		  WRITE(7,57) K,L
57		  FORMAT(/,' FOR GROUP ',I2, X, ' AND PARTITION',I2,X,
     :		  ' THE MEAN IS')
		  WRITE(7,58)(EMU(K,L,J),J=1,IPC(L))
58		  FORMAT(10F12.6)
		ELSE
		  WRITE(7,57) K,L
		  WRITE(7,58)((EMUL(K,L,J,M),M=1,IM(L)),J=1,IPC(L))
		  DO J=ISV(L),IEV(L)	!62
		   IF (IVARTYPE(J).EQ.3) THEN
		    WRITE(7,54)
		    WRITE(7,55)(THETA(K,J,M),M=1,NCAT(J))
		   END IF
  		  END DO	!62

	   END IF
  	END DO	!56
  	END DO	!56
	DO K=1,NG	!59
	DO L=1,NPAR	!59
	     IF (IPARTYPE(L).NE.1) THEN
		WRITE(7,60) L,K
60		FORMAT(/,2X,'VARIANCE FOR PARTITION',I2,X,' AND GROUP',I2)
		DO I=1,IPC(L)	!63
		   WRITE(7,61)(VARIX(K,L,I,J),J=1,IPC(L))
61		   FORMAT(10F13.6)
  		END DO	!63

	     END IF
  	END DO	!59
  	END DO	!59

*	Print out the types of models for each partition

	DO L=1,NPAR	!74
	   IF (IPARTYPE(L).EQ.1) THEN
		WRITE(7,75) L,IP(L),IPARTYPE(L)
75		FORMAT(' PARTITION',I3,' HAS',I2,' VARIABLES',/,' ITYPE',
     :   ' IS',I2,' HENCE A CATEGORICAL MODEL FOR THIS PARTITION')
	   ELSE IF (IPARTYPE(L).EQ.2) THEN
		WRITE(7,76) L,IP(L),IPARTYPE(L)
76		FORMAT(' PARTITION',I3,' HAS',I2,' VARIABLES',/,' ITYPE',
     :   ' IS',I2,' HENCE A MVN MODEL FOR THIS PARTITION')
	   ELSE
		WRITE(7,77) L,IP(L),IPARTYPE(L)
77		FORMAT(' PARTITION',I3,' HAS',I2,'VARIABLES',/,' ITYPE',
     :   ' IS',I2,' HENCE A LOCATION MODEL FOR THIS PARTITION')
	END IF
  	END DO	!74


*	E step of the EM algorithm.
*	Estimate the complete data sufficient statistics given the
*	data, & current values of means,variances and mixing proportions.

	ICNT=1

*	Form the augmented covariance matrix

99999	DO K=1,NG	!40
	DO L=1,NPAR	!40
	    IF (IPARTYPE(L).EQ.2) THEN
		AUGCOV(K,L,0,0)=-1.
		DO J=1,IPC(L)	!41
		    AUGCOV(K,L,J,0)=EMU(K,L,J)
		    AUGCOV(K,L,0,J)=AUGCOV(K,L,J,0)
		    DO J1=1,J	!41
		    AUGCOV(K,L,J,J1)=VARIX(K,L,J,J1)
		    AUGCOV(K,L,J1,J)=AUGCOV(K,L,J,J1)
  		END DO	!41
  		END DO	!41
	    END IF
  	END DO	!40
  	END DO	!40


	ALL=0.0
	ALIM=1E-7
	DO II = 1,NOBS	!20
	   SDENS = 0.0
	   DO K=1,NG	!21
	      PROB(K,II)=1.0
	      APRODENS(K)=1.0

*	evaluate the discrete variables contribution to the densities

	      DO L=1,NPAR	!37
		IF (IPARTYPE(L).EQ.1) THEN
		   DO J=ISV(L),IEV(L)	!22
		   DO M=1,NCAT(J)	!22
			IF (IIX(II,J,M).EQ.1) THEN
			   PROB(K,II)=THETA(K,J,M)*PROB(K,II)
			END IF
  		   END DO	!22
  		   END DO	!22
		ELSE IF (IPARTYPE(L).EQ.3) THEN
		   DO J=ISV(L),IEV(L)	!35
			IF (IVARTYPE(J).EQ.3) THEN
			   DO M=1,NCAT(J)	!36
			     IF (IIX(II,J,M).EQ.1) THEN
			   	PROB(K,II)=THETA(K,J,M)*PROB(K,II)
*			ELSE IF (IIX(II,J,M).EQ.9) THEN
*			   IF (THETA(K,J,M).GT.ALIM
*    :			     PROB(K,II)=(THETA(K,J,M)**THETA(K,J,M))
*     :						*PROB(K,II)
			     END IF
  			   END DO	!36

			END IF
  		   END DO	!35

		END IF
  	      END DO	!37


*	*******HAVEN'T DONE THE MISSING PART FOR THE LOCATION MODEL
*	*******CHECK CATEGORICAL CONTRIBNTO LIKELIHOOD FOR MISSING***

*	evaluate the continuous variables contribution to the densities

	      DO L=1,NPAR	!23
		IF (IPARTYPE(L).NE.1) THEN
		  DENS(K,L)=0.0
		  I1=0

*	check whether variables are missing in a partition, and create
*	ISWPCOL(j). The augmented covariance matrix is swept on the
*	columns specified in ISWPCOL.

		IF (IPARTYPE(L).EQ.2) THEN
		   ISUMISS=0
		   J1=0
		   DO J=ISV(L),IEV(L)	!42
			J1=J1+1
			ISUMISS=ISUMISS + IMISSCODE(II,J)
			ISWPCOL(J1) = IMISSCODE(II,J)
  		   END DO	!42


*	If some variables are missing in the partition, sweep and then
*	do the regression

		IF ((ISUMISS.GT.0).AND.(ISUMISS.LE.IPC(L))) THEN
     		  CALL SWEEP(AUGCOV,ISWPCOL,ISUMISS,K,L,IPC(L),EN)
		  J3 = 0
		  DO J = ISV(L),IEV(L)	!43
		     J3=J3+1
		     IF (IMISSCODE(II,J).EQ.1) THEN
			ESTX(K,II,J)=EN(K,L,0,J3)
			J2=0
			DO I=ISV(L),IEV(L)	!44
			   J2=J2+1
			   IF (IMISSCODE(II,I).EQ.0) ESTX(K,II,J) =
     :				ESTX(K,II,J) + EN(K,L,J3,J2)*X(II,I)
  			END DO	!44

		     END IF
  		  END DO	!43

		END IF
		END IF

*	Look at the contribution of the observed values to the likelihood
*	First, the inverse and the determinant of the covariance matrix
*	for the observed values of observation I is found by using NAG.

	IF (IPARTYPE(L).EQ.2) THEN
	   I1=0
	   I3=0
	   DO I=ISV(L),IEV(L)	!45
		I3=I3+1
		I2=0
		IF (IMISSCODE(II,I).EQ.0) THEN
		   I1=I1+1
		   I4=0
		   DO J=ISV(L),IEV(L)	!46
			I4=I4+1
			IF(IMISSCODE(II,J).EQ.0) THEN
			   I2=I2+1
			   VAROBS(K,L,I1,I2)=VARIX(K,L,I3,I4)
			END IF
  		   END DO	!46

		END IF
  	   END DO	!45

	NDIM=IPC(L)-ISUMISS
	IF (ISUMISS.EQ.0) THEN
	   CALL DETINV(VAROBS,K,L,IPC(L),ADET,VARIN)
	ELSE IF (ISUMISS.LT.IPC(L)) THEN
	   CALL DETINV(VAROBS,K,L,NDIM,ADET,VARIN)
	END IF
	END IF

*	(i) evaluate the MVN contribution
		  IF(IPARTYPE(L).EQ.2) THEN
			I1=0
			I2=0
			DO I=ISV(L),IEV(L)	!24
			    J1=0
			    I1=I1+1
			    IF (IMISSCODE(II,I).EQ.0) THEN
				I2=I2+1
				J2=0
			  	DO J=ISV(L),IEV(L)	!47
				   J1=J1+1
				   IF (IMISSCODE(II,J).EQ.0) THEN
					J2=J2+1
	           			DENS(K,L) = DENS(K,L)
     :    + (X(II,I)-EMU(K,L,I1))*VARIN(K,L,I2,J2)*(X(II,J)-EMU(K,L,J1))
				   END IF
  				END DO	!47

			    END IF
  			  END DO	!24

*	(ii) evaluate the continuous location variables contribution
		  ELSE IF(IPARTYPE(L).EQ.3) THEN
			DO I=ISV(L),IEV(L)	!25
			 IF(IVARTYPE(I).EQ.3) M=IX(II,I)
  			END DO	!25

			DO I=ISV(L),IEV(L)	!26
			 IF (IVARTYPE(I).EQ.4) THEN
			   J1=0
			   I1=I1+1
			   DO J=ISV(L),IEV(L)	!30
			    IF (IVARTYPE(J).EQ.4) THEN
			     J1=J1+1
			     DENS(K,L)=DENS(K,L) + (X(II,I)-EMUL(K,L,I1,M))
     :  *VARIN(K,L,I1,J1)*(X(II,J)-EMUL(K,L,J1,M))
			    END IF
  			   END DO	!30

			 END IF
  			END DO	!26

		  END IF
	IF (ISUMISS.LT.IPC(L)) THEN
		  DENS(K,L) = DEXP(-0.5*DENS(K,L))
		  A=0.5*FLOAT(IPC(L)-ISUMISS)
		  DENS(K,L)=DENS(K,L)/((2.0*PIE)**(A)*DSQRT(ADET(K,L)))
		  APRODENS(K)=APRODENS(K)*DENS(K,L)
	END IF
		END IF
  	     END DO	!23

	     APRODENS(K)=APRODENS(K)*PROB(K,II)*PI(K)
	     SDENS=SDENS+APRODENS(K)
  	   END DO	!21

	   IF (SDENS.NE.0.0) THEN
		DO K=1,NG	!27
		    Z(II,K)=APRODENS(K)/SDENS
  		END DO	!27

		ALL=ALL+DLOG(SDENS)
	   ELSE
		DO K=1,NG	!28
		   Z(II,K)=0.01
  		END DO	!28

		   WRITE(7,29)II
29		   FORMAT(//' SUM OF DENSITY FUNCTIONS IS ZERO',
     :		   ' FOR OBSERVATION',I4)
	   END IF
  	END DO	!20


*	Check on convergence - look at the likelihood function
*	If the absolute value of the difference in 2 likelihoods is
*	less than a tolerance value the estimates are written out.
*	A check is made on the number of iterations (200 max).
*	statement 100 - the M step
*	statement 500 - print out the current estimates (algorithm
*	has converged)
*	statement 501 - algorithm hasn't converged, current estimates
*	are printed out.

	CLOGLI(ICNT) = ALL
	WRITE(7,888),ICNT,CLOGLI(ICNT)
888	FORMAT(1X,'FOR LOOP',I5,'LOGLIKELIHOOD IS',F16.8)
	IF (ICNT.LE.10) GO TO 100
	C=1.D-10
	TOL=ABS(CLOGLI(ICNT) - CLOGLI(ICNT-10))
	IF(TOL.LE.C) GO TO 500
	IF(ICNT.GE.iter) GO TO 501

*	M step of the EM algorithm
*	Calculate an updated estimate of the parameters
*	 means and covariances.

*	(i) the mixing proportions,

100 	DO K=1,NG	!101
	 ZSUM(K)=0.0
	   DO II=1,NOBS	!102
	  	ZSUM(K)=ZSUM(K) + Z(II,K)
   	   END DO	!102

	 PI(K) = ZSUM(K)/NOBS
     	END DO	!101


*	(ii) the conditional probabilities

	DO K=1,NG	!103
	DO L=1,NPAR	!103
		IF(IPARTYPE(L).EQ.1) THEN
		  DO J=ISV(L),IEV(L)	!104
		  DO M=1,NCAT(J)	!104
			  THETA(K,J,M) = 0.0
			  DO II=1,NOBS	!105
			     IF(IMISSCODE(II,J).EQ.1) THEN
				THETA(K,J,M) = THETA(K,J,M) +
     :				Z(II,K)*THETA2(K,J,M)
			     ELSE IF(IMISSCODE(II,J).EQ.0) THEN
				THETA(K,J,M) = THETA(K,J,M) +
     :				 Z(II,K)*IIX(II,J,M)
			     END IF
   			  END DO	!105

			  THETA(K,J,M)=THETA(K,J,M)/ZSUM(K)
   		  END DO	!104
   		  END DO	!104
		ELSE IF(IPARTYPE(L).EQ.2) THEN

*	(iii) the means (EMU(K,L,J))

		  J1=0
		  DO J=ISV(L),IEV(L)	!106
		   J1=J1+1
		   XSUM(K,J1)=0.0
		   DO II=1,NOBS	!107
		      IF(IMISSCODE(II,J).EQ.1) X(II,J)=ESTX(K,II,J)
		      XSUM(K,J1)=XSUM(K,J1) + X(II,J)*Z(II,K)
   		   END DO	!107

		   EMU(K,L,J1)=XSUM(K,J1)/ZSUM(K)
   		  END DO	!106

*	(iv) the location model parameters
*	   (i)the discrete variable

		ELSE IF(IPARTYPE(L).EQ.3) THEN
		  DO J=ISV(L),IEV(L)	!108
		    IF(IVARTYPE(J).EQ.3) THEN
			DO M=1,NCAT(J)	!109
			  THETA(K,J,M)=0.0
			  DO II=1,NOBS	!110
			     IF(IX(II,J).EQ.M) THETA(K,J,M)
     :					= THETA(K,J,M) + Z(II,K)
   			  END DO	!110

			  THETA(K,J,M)=THETA(K,J,M)/ZSUM(K)
   		        END DO	!109

		    END IF
   		  END DO	!108


*		(ii) the continuous variables the means EMUL(K,L,J,M)
		DO J=ISV(L),IEV(L)	!111
		  IF (IVARTYPE(J).EQ.3) THEN
		     DO M=1,IM(L)	!112
			J1=0
			DO JJ=ISV(L),IEV(L)	!112
			  IF (IVARTYPE(JJ).EQ.4) THEN
			    J1=J1+1
			    XSUM2(K,J1,M)=0.0
		            ZM(K,M)=0.0
			    DO II=1,NOBS	!113
			      IF(IX(II,J).EQ.M) THEN
				XSUM2(K,J1,M)=XSUM2(K,J1,M) +
     :					     X(II,JJ)*Z(II,K)
				ZM(K,M)=ZM(K,M)+Z(II,K)
			      END IF
   			     END DO	!113

			  EMUL(K,L,J1,M)=XSUM2(K,J1,M)/ZM(K,M)
			  END IF
   			END DO	!112
   			END DO	!112
		    END IF
   		END DO	!111

		WRITE(7,601)K,((EMUL(K,L,J1,M),M=1,IM(L)),J1=1,IPC(L))
601	FORMAT(/2X,'FOR GROUP',I2,X,'EMUL',10F10.4)
		END IF
   	END DO	!103
   	END DO	!103


*	Calculate updated estimates of the variances
*	(i) the multivariate normal data

	DO K=1,NG	!115
		DO L=1,NPAR	!116
		 IF (IPARTYPE(L).NE.1) THEN
		   DO J=1,IPC(L)	!117
	 	   DO I=1,J	!117
		        VAR(K,L,I,J)=0.0
   		   END DO	!117
   		   END DO	!117
		 END IF
   		END DO	!116

		DO II=1,NOBS	!118
		DO L=1,NPAR	!118
		      IF (IPARTYPE(L).EQ.2) THEN
			J1=0
			DO J=ISV(L),IEV(L)	!119
		           IF(IMISSCODE(II,J).EQ.1) X(II,J)=ESTX(K,II,J)
			   I1=0
			   J1=J1+1
			   DO I=ISV(L),J	!119
		      	      IF(IMISSCODE(II,I).EQ.1)X(II,I)=ESTX(K,II,I)
			      I1=I1+1
			      IF ((IMISSCODE(II,I).EQ.1) .AND.
     :				(IMISSCODE(II,J).EQ.1)) THEN
				VAR(K,L,I1,J1)=VAR(K,L,I1,J1) +
     :     ((X(II,J) - EMU(K,L,J1))*(X(II,I)-EMU(K,L,I1))
     :                + EN(K,L,I1,J1))*Z(II,K)
			      ELSE
			        VAR(K,L,I1,J1) = VAR(K,L,I1,J1) + (X(II,J)
     :			        -EMU(K,L,J1))*(X(II,I)-EMU(K,L,I1))*Z(II,K)
			      END IF
   			END DO	!119
   			END DO	!119
		   END IF
   		   END DO	!118
   		   END DO	!118

*	(ii) the continuous location data

		   DO L=1,NPAR	!120
		      IF (IPARTYPE(L).EQ.3) THEN
			DO J=ISV(L),IEV(L)	!121
			 IF (IVARTYPE(J).EQ.3) THEN
			  DO II=1,NOBS	!122
			  DO M=1,IM(L)	!122
			    IF(IX(II,J).EQ.M) THEN
			     J1=0
			     DO JJ=ISV(L),IEV(L)	!123
				IF (IVARTYPE(JJ).EQ.4) THEN
			         I1=0
			         J1=J1+1
			          DO I=ISV(L),J	!124
			           IF (IVARTYPE(I).EQ.4) THEN
			            I1=I1+1
			      VAR(K,L,I1,J1) = VAR(K,L,I1,J1)+(X(II,JJ)
     :			      -EMUL(K,L,J1,M))*(X(II,I)-EMUL(K,L,I1,M))
     :					*Z(II,K)
			    	   END IF
   				  END DO	!124

				END IF
   			     END DO	!123

			    END IF
   			  END DO	!122
   			  END DO	!122
			 END IF
   			END DO	!121

		      END IF
   		   END DO	!120

		DO L=1,NPAR	!126
		 IF (IPARTYPE(L).NE.1) THEN
		   DO J=1,IPC(L)	!125
		   DO I=1,J	!125
			VAR(K,L,I,J)=VAR(K,L,I,J)/ZSUM(K)
			VAR(K,L,J,I)=VAR(K,L,I,J)
   		   END DO	!125
   		   END DO	!125
		 END IF
   		END DO	!126

   	END DO	!115


*	Make a copy of the covariance matrix before we use NAG

	DO K=1,NG	!130
	DO L=1,NPAR	!130
	    IF(IPARTYPE(L).NE.1) THEN
		DO J=1,IPC(L)	!131
		DO I=1,IPC(L)	!131
		   VARIX(K,L,I,J)=VAR(K,L,I,J)
   	        END DO	!131
   	        END DO	!131
	    END IF
   	END DO	!130
   	END DO	!130

*	Make a copy of the conditional probabilities for use with
*	missing categorical data

	DO K=1,NG	!132
	DO J=1,NVAR	!132
	   IF (IVARTYPE(J).EQ.1) THEN
		DO M=1,NCAT(J)	!133
		   THETA2(K,J,M)=THETA(K,J,M)
   		END DO	!133

	   END IF
   	END DO	!132
   	END DO	!132

	ICNT=ICNT+1

*	send back to the E step
	GO TO 99999

*	Write out the current estimates of the parameters.

*	(1) If the algorithm has not converged:-

501	WRITE(7,502)
502	FORMAT(/'----------------------------------------------------')
	WRITE(7,503) iter
503	FORMAT(/' THE EM ALGORITHM HAS NOT CONVERGED AFTER ',I3,
     :  /,' ITERATIONS BUT THE CURRENT ESTIMATES OF THE PARAMETERS ',
     :	/,' WILL BE PRINTED OUT.')
	WRITE(7,502)

*	The parameters are to be written to 'EMPARAMEST.DAT' to be
*	used as input for the PROGRAM MULTIMIX. ISPEC is set to 1.

	OPEN(11, FILE='EMPARAMEST.OUT',STATUS='NEW')
	ISPEC=1
	WRITE(11,504) NG, NOBS, NVAR, NPAR, ISPEC
504	FORMAT(X,5I6)
	WRITE(11,505) (JP(J),J=1,NVAR)
	WRITE(11,505) (IP(L),L=1,NPAR)
	WRITE(11,505) (IPC(L),L=1,NPAR)
	WRITE(11,505) (ISV(L),L=1,NPAR)
	WRITE(11,505) (IEV(L),L=1,NPAR)
	WRITE(11,505) (IPARTYPE(L),L=1,NPAR)
	WRITE(11,505) (IVARTYPE(J),J=1,NVAR)
	WRITE(11,505) (NCAT(J),J=1,NVAR)
505	FORMAT(10I4)
	WRITE(11,506) (PI(K),K=1,NG)
506	FORMAT(10F10.6)
	DO K=1,NG	!507
	DO L=1,NPAR	!507
	   IF (IPARTYPE(L).EQ.1) THEN
	      DO J=ISV(L),IEV(L)	!508
		WRITE(11,509)(THETA(K,J,M),M=1,NCAT(J))
509		FORMAT(10F10.6)
   	      END DO	!508

	   ELSE IF(IPARTYPE(L).EQ.2) THEN
	      WRITE(11,510) (EMU(K,L,J),J=1,IPC(L))
510	      FORMAT(10F13.6)
	   ELSE
	      DO J=ISV(L),IEV(L)	!511
		IF (IVARTYPE(J).EQ.3) THEN
		  WRITE(11,509)(THETA(K,J,M),M=1,NCAT(J))
		END IF
   	      END DO	!511

	      DO J=1,IPC(L)	!512
		WRITE(11,510) (EMUL(K,L,J,M),M=1,IM(L))
   	      END DO	!512

	   END IF
   	END DO	!507
   	END DO	!507
	DO K=1,NG	!513
	DO L=1,NPAR	!513
	   IF (IPARTYPE(L).NE.1) THEN
	      DO I=1,IPC(L)	!514
		 WRITE(11,515)(VARIX(K,L,I,J),J=1,IPC(L))
515		 FORMAT(10F13.6)
   	      END DO	!514

	   END IF
   	END DO	!513
   	END DO	!513

*	(2) the current estimates of the parameters are printed out
*	Estimates of the proportions in each group and loglikelihood

500	WRITE(7,888),ICNT,CLOGLI(ICNT)
	DO K=1,NG	!540
	   WRITE(7,541) K,PI(K)
541	   FORMAT(/' THE ESTIMATE OF THE MIXING PROPORTION IN GROUP ',I3,
     :	   'IS ', F10.8)
   	END DO	!540


*	estimates of the probabilities for each group

	DO L=1,NPAR	!525
	  WRITE(7,502)
	  WRITE(7,547)L
547	FORMAT(/,' THE CURRENT ESTIMATES FOR PARTITION',I3)
	DO K=1,NG	!525
	    WRITE(7,524)K
524	    FORMAT(/,' Group:',I3,/,' ---------')
	   IF (IPARTYPE(L).EQ.1) THEN
		DO J=ISV(L),IEV(L)	!528
		  WRITE(7,529)J
529		  FORMAT(' For variable ',I3,X,'the level ',
     :	'probabilities are')
		  WRITE(7,530) (THETA(K,J,M),M=1,NCAT(J))
530		  FORMAT(10F10.6)
   	        END DO	!528

	   ELSE IF(IPARTYPE(L).EQ.2) THEN

*	Estimate of the means for each group.

		WRITE(7,531) (EMU(K,L,J),J=1,IPC(L))
531		FORMAT(' The mean for this partition is ',/,10F13.6)

*	estimates of the location model parameters

	ELSE IF(IPARTYPE(L).EQ.3) THEN
		DO J=ISV(L),IEV(L)	!533
		  IF (IVARTYPE(J).EQ.3) THEN
		   WRITE(7,534)J,(THETA(K,J,M),M=1,NCAT(J))
534		   FORMAT(' For variable',I3,X,'THETA(K,J,M)'
     :					' is',10F10.6)
		  END IF
   		END DO	!533

		J1=0
		DO J=ISV(L),IEV(L)	!535
		 IF(IVARTYPE(J).EQ.4) THEN
		   J1=J1+1
		   WRITE(7,536) (EMUL(K,L,J1,M),M=1,IM(L))
536		   FORMAT(' The mean for the continuous location'
     :				' variables is',/,10F13.6)
		 END IF
   		END DO	!535

	END IF
   	END DO	!525
   	END DO	!525

*	Estimates of the variances for each group.

	DO L=1,NPAR	!523
	   IF(IPARTYPE(L).NE.1) THEN
	        WRITE(7,502)
		WRITE(7,522)L
522		FORMAT(/,' THE CURRENT ESTIMATE OF THE COVARIANCE MATRIX FOR',
     :  ' PARTITION ',I3)
	DO K=1,NG	!543
	   WRITE(7,526)K
526	   FORMAT(/,' Group',I3,':',/,' --------')
	   DO I=1,IPC(L)	!543
		   WRITE(7,527)(VAR(K,L,I,J),J=1,IPC(L))
527		   FORMAT(10F15.6)
   	   END DO	!543
   	   END DO	!543
	      END IF
   	END DO	!523

	WRITE(7,502)

*	Determine the assignment of the observations to groups

	DO I=1,NOBS	!516
		IMAX = 1
		   DO K=2,NG	!517
		      IF(Z(I,K).GT.Z(I,IMAX)) THEN
		         IMAX=K
		      END IF
   		   END DO	!517

		IGP(I)=IMAX
   	END DO	!516

	WRITE(7,542)
542	FORMAT(/,X,'THE ASSIGNMENT OF OBSERVATIONS TO GROUPS',/)
	WRITE(7,518) (IGP(I),I=1,NOBS)
518	FORMAT(10I3)
	DO K=1,NG	!537
		NUM(K)=0
   	END DO	!537

	DO I=1,NOBS	!538
	DO K=1,NG	!538
		IF(IGP(I).EQ.K) NUM(K)=NUM(K)+1
   	END DO	!538
   	END DO	!538
	WRITE(7,539)(NUM(K),K=1,NG)
539	FORMAT(1X,/,' TOTAL NUMBERS IN EACH GROUP',/,10I5)

*	The estimates of the missing data given the group assignment

	WRITE(7,502)
	WRITE(7,544)
544	FORMAT(/,' THE ESTIMATES OF THE MISSING DATA VALUES GIVEN THE ',
     :  'GROUP ASSIGNMENT',//,
     :	5X,'GROUP',5X,'OBSERVATION',4X,'VARIABLE',3X,'ESTIMATED VALUE')
	DO II=1,NOBS	!545
	DO J=1,NVAR	!545
	DO K=1,NG	!545
	    IF (IMISSCODE(II,J).EQ.1) THEN
	      IF (IGP(II).EQ.K) THEN
		IF (IVARTYPE(J).EQ.2) THEN
		   WRITE(7,546)K,II,J,ESTX(K,II,J)
546		   FORMAT(7X,I2,8X,I5,10X,I2,6X,F16.8)
		END IF
	      END IF
	    END IF
   	END DO	!545
   	END DO	!545
   	END DO	!545
	WRITE(7,502)

*	have a look at the Zij,s, and write out the assigned groups and
*	the Zij's to GROUPS.OUT

	WRITE(7,521)
521	FORMAT(//' THE ESTIMATES OF THE POSTERIOR PROBALITIES')
		DO I=1,NOBS	!519
		   WRITE(7,520)I,(Z(I,K),K=1,NG)
520		   FORMAT(' OBSERVATION',I4,2X,10F10.6)
		   WRITE(12,532) IGP(I), (Z(I,K),K=1,NG)
532		   FORMAT(X,I2,X,10F9.6)
   		END DO	!519

	END

*	This subroutine calculates the inverse and the determinant of
*	the covariance matrix using NAG routines F03ABF and F01ABF.

	SUBROUTINE DETINV(VAROBS,IK,IL,IDIM,ADET,VARIN)
	IMPLICIT DOUBLE PRECISION(A-H,O-Z)
	PARAMETER (ik6=6,ip15=15,IIP=ip15+1)
	DIMENSION VAROBS(ik6,ip15,ip15,ip15),VARIN(ik6,ip15,ip15,ip15),
     :  ADET(ik6,ip15)
	INTEGER IA,IFAIL,IB
	REAL*8 TEMP(IIP,IIP),DET,WKSPCE(IIP),ZI(IIP),B(IIP,IIP)
	IA=IIP
	IB=IIP
	IFAIL = 0
	    	N=IDIM
  	    	DO I=1,IDIM	!202
	    	DO J=1,IDIM	!202
		   TEMP(I,J) = VAROBS(IK,IL,I,J)
   	    	END DO	!202
   	    	END DO	!202
		CALL F03ABF(TEMP,IA,N,DET,WKSPCE,IFAIL)
		CONTINUE
		IF (IFAIL.EQ.0) THEN
		   ADET(IK,IL)=DET
		ELSE
		   WRITE(7,204)IFAIL,IK,IL
204		   FORMAT(//,' TO CALCULATE THE DETERMINANT IFAIL = ',I2,
     :		   ' FOR GROUP',I3,' PARTITION',I3)
		END IF
		IFAIL=0
		CALL F01ABF(TEMP,IA,N,B,IB,ZI,IFAIL)
		CONTINUE
		IF(IFAIL.EQ.0) THEN
		   DO I = 1,IDIM	!213
		   DO J = 1,IDIM	!213
		      VARIN(IK,IL,I,J)=B(I,J)
   		   END DO	!213
   		   END DO	!213
		ELSE
		   WRITE(7,207)IFAIL,IK,IL
207		   FORMAT(//,' TO CALCULATE THE INVERSE IFAIL = ',I2,
     :		   ' FOR GROUP', I2,' PARTITION', I3)
		END IF
	DO J=1,IDIM	!199
	     DO I=1,J	!199
	        VARIN(IK,IL,I,J)=VARIN(IK,IL,J,I)
   	END DO	!199
   	END DO	!199
	RETURN
	END

*	This subroutine sweeps the augmented covariance matrix on the
*	observed values. Columns are specified in ISWPCOL
*	EN = SWP [i1,i2,...,it] EM where i1,i2,..,it are the columns
*	specified in ISWPCOL

	SUBROUTINE SWEEP(AUG,ISP,ISUM,IK,IL,IC,EN)
	IMPLICIT DOUBLE PRECISION(A-H,O-Z)
	PARAMETER(ik6=6,ip15=15)
	DIMENSION EM(ik6,ip15,0:ip15,0:ip15),EN(ik6,ip15,0:ip15,0:ip15),
     :		ISP(ip15),AUG(ik6,ip15,0:ip15,0:ip15)

*	First we make a copy of the augmented covariance matrix

	DO I=0,IC	!315
	DO J=0,IC	!315
	   EM(IK,IL,I,J)=AUG(IK,IL,I,J)
   	END DO	!315
   	END DO	!315

	IF (ISUM.LT.IC) THEN
	  DO J1=1,IC	!303
	    IF (ISP(J1).EQ.0) THEN
		I=J1

*	Calculation of the diagonal elements

		EN(IK,IL,I,I) = -1/EM(IK,IL,I,I)

*	Remainder of the iTH row and column

		DO J=0,IC	!300
		  IF(J.NE.I) THEN
			EN(IK,IL,I,J)=EM(IK,IL,I,J)/EM(IK,IL,I,I)
			EN(IK,IL,J,I)=EN(IK,IL,I,J)
		  END IF
   		END DO	!300


*	The remaining elements of N

		DO J=0,IC	!301
	    	   IF (J.NE.I) THEN
			DO K=0,J	!302
		   	   IF(K.NE.I) THEN
				EN(IK,IL,J,K) = EM(IK,IL,J,K)
     :  				- EM(IK,IL,J,I) * EN(IK,IL,I,K)
				EN(IK,IL,K,J)=EN(IK,IL,J,K)
		   	   END IF
   			END DO	!302

	    	   END IF
   		END DO	!301

		DO JJ=0,IC	!305
		DO II=0,IC	!305
		    EM(IK,IL,II,JJ)=EN(IK,IL,II,JJ)
   		END DO	!305
   		END DO	!305
	    END IF
   	  END DO	!303

	ELSE IF(ISUM.EQ.IC) THEN
	  DO J=0,IC	!304
	  DO J2=0,J	!304
	      EN(IK,IL,J,J2)=EM(IK,IL,J,J2)
	      EN(IK,IL,J2,J)=EN(IK,IL,J,J2)
   	  END DO	!304
   	  END DO	!304
	END IF
	RETURN
	END
