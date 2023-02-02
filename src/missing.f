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
	  DO 12 K=1,NG
	  DO 12 L=1,NPAR
		IF (IPARTYPE(L).EQ.1) THEN
			DO  13 J=ISV(L),IEV(L)
			  READ(9,*) (THETA(K,J,M),M=1,NCAT(J))
13			CONTINUE
		ELSE IF(IPARTYPE(L).EQ.2) THEN
			READ(9,*)(EMU(K,L,J),J=1,IPC(L))
		ELSE IF(IPARTYPE(L).EQ.3) THEN
			DO 14 J=ISV(L),IEV(L)
			  IF(IVARTYPE(J).EQ.3) THEN 
				READ(9,*)(THETA(K,J,M),M=1,NCAT(J))
				IM(L)=NCAT(J)
			   END IF
14			CONTINUE
			J1=0
			DO 15 J=ISV(L),IEV(L)
			  IF(IVARTYPE(J).EQ.4) THEN
			   J1=J1+1
			   READ(9,*) 
     :			    (EMUL(K,L,J1,M),M=1,IM(L))
			  END IF
15			CONTINUE
		END IF
12	CONTINUE

	DO 16 K=1,NG
	DO 16 L=1,NPAR
		IF(IPARTYPE(L).NE.1) 
     :		READ(9,*)((VARIX(K,L,I,J),J=1,IPC(L)),I=1,IPC(L))
16  	CONTINUE

	ELSE
	  READ(9,*) (IGRP(I),I=1,NOBS)
	  DO 18 I=1,NOBS
	  DO 18 K=1,NG
	      Z(I,K)=0.0
18	  CONTINUE
	  DO 31 I=1,NOBS
	    IK=IGRP(I)
	    Z(I,IK)=1.0
31	  CONTINUE
	  DO 33 L=1,NPAR
	     IF(IPARTYPE(L).EQ.3) THEN
		DO 32 J=ISV(L),IEV(L)
		  IF (IVARTYPE(J).EQ.3) IM(L)=NCAT(J)
32		CONTINUE
	    END IF
33	  CONTINUE
	END IF

*	read in the observations 

	DO 17 I=1,NOBS
		READ(8,* ) (X(I,JP(J)),J=1,NVAR)
17	CONTINUE

*	Check for the presence of missing values, and create 
*	IMISSCODE(I,J), where IMISSCODE(I,J)=1 if observation i variable j is 
*	missing, and IMISSCODE(I,J)=0 otherwise

	DO 34 I=1,NOBS
	DO 34 J=1,NVAR
		IVAL=NINT(X(I,J))
		IF (IVAL.EQ.IMISS) THEN
			IMISSCODE(I,J)=1
		ELSE
			IMISSCODE(I,J)=0
		END IF
34	CONTINUE

*	Check the categorical variable to see if 0 is a category.
*	Add 1 to the variables if it is so that Xij=0 is taken as 
*	level 1, Xij=1 is level 2 etc.

	DO 80 J=1,NVAR
	    IF((IVARTYPE(J).EQ.1).OR.(IVARTYPE(J).EQ.3)) THEN
		IMIN=NINT(X(1,J))
		IF (IMIN.NE.0) THEN
		   DO 81 I=2,NOBS
			INEXT=NINT(X(I,J))
			IF(INEXT.EQ.0) IMIN=INEXT
81		   CONTINUE
		END IF
		IF (IMIN.EQ.0) THEN
		    DO 82 I=1,NOBS
			IF (IMISSCODE(I,J).EQ.0) THEN
		 	   IX(I,J)=NINT(X(I,J)) + 1
			ELSE
			   IX(I,J)=NINT(X(I,J))
			END IF
82		    CONTINUE
		ELSE
		    DO 83 I=1,NOBS
			IX(I,J)=NINT(X(I,J))
83		    CONTINUE
		END IF
	   END IF
80	CONTINUE

*	For the categorical variables, create the matrix IIX(I,J,M) where 
*	x(i,j,m)=1 if obs. i, variable j is at level m,
*	x(i,j,m)=0 if not at level m.

	DO 84 I=1,NOBS
	DO 84 J=1,NVAR
	    IF((IVARTYPE(J).EQ.1).OR.(IVARTYPE(J).EQ.3)) THEN
		DO 85 M=1,NCAT(J)
		     IF (IX(I,J).EQ.M) THEN
			IIX(I,J,M)=1
		     ELSE IF (IX(I,J).EQ.IMISS) THEN
			IIX(I,J,M)=9
		     ELSE
			IIX(I,J,M)=0
		     END IF
85		CONTINUE
	     END IF
84	CONTINUE

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
	DO 56 K=1,NG
	DO 56 L=1,NPAR
		IF(IPARTYPE(L).EQ.1) THEN
		  WRITE(7,54) K
54		  FORMAT(/,2X,' THETA(K,J,M) FOR GROUP ',I2)
		  DO 64 J=ISV(L),IEV(L)
		      WRITE(7,55)(THETA(K,J,M),M=1,NCAT(J))
55		      FORMAT(10F8.4)
64		  CONTINUE
		ELSE IF(IPARTYPE(L).EQ.2) THEN
		  WRITE(7,57) K,L
57		  FORMAT(/,' FOR GROUP ',I2, X, ' AND PARTITION',I2,X, 
     :		  ' THE MEAN IS')
		  WRITE(7,58)(EMU(K,L,J),J=1,IPC(L))
58		  FORMAT(10F12.6)
		ELSE
		  WRITE(7,57) K,L
		  WRITE(7,58)((EMUL(K,L,J,M),M=1,IM(L)),J=1,IPC(L))
		  DO 62 J=ISV(L),IEV(L)
		   IF (IVARTYPE(J).EQ.3) THEN
		    WRITE(7,54)
		    WRITE(7,55)(THETA(K,J,M),M=1,NCAT(J))
		   END IF
62		  CONTINUE
	   END IF
56	CONTINUE
	DO 59 K=1,NG
	DO 59 L=1,NPAR
	     IF (IPARTYPE(L).NE.1) THEN
		WRITE(7,60) L,K
60		FORMAT(/,2X,'VARIANCE FOR PARTITION',I2,X,' AND GROUP',I2)
		DO 63 I=1,IPC(L)
		   WRITE(7,61)(VARIX(K,L,I,J),J=1,IPC(L))
61		   FORMAT(10F13.6)
63		CONTINUE
	     END IF
59	CONTINUE
 
*	Print out the types of models for each partition

	DO 74 L=1,NPAR
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
74	CONTINUE

*	E step of the EM algorithm.
*	Estimate the complete data sufficient statistics given the
*	data, & current values of means,variances and mixing proportions.

	ICNT=1

*	Form the augmented covariance matrix

99999	DO 40 K=1,NG
	DO 40 L=1,NPAR
	    IF (IPARTYPE(L).EQ.2) THEN
		AUGCOV(K,L,0,0)=-1.
		DO 41 J=1,IPC(L)
		    AUGCOV(K,L,J,0)=EMU(K,L,J)
		    AUGCOV(K,L,0,J)=AUGCOV(K,L,J,0)
		    DO 41 J1=1,J
		    AUGCOV(K,L,J,J1)=VARIX(K,L,J,J1)
		    AUGCOV(K,L,J1,J)=AUGCOV(K,L,J,J1)
41		CONTINUE
	    END IF	
40	CONTINUE


	ALL=0.0
	ALIM=1E-7
	DO 20 II = 1,NOBS
	   SDENS = 0.0
	   DO 21 K=1,NG
	      PROB(K,II)=1.0
	      APRODENS(K)=1.0

*	evaluate the discrete variables contribution to the densities

	      DO 37 L=1,NPAR
		IF (IPARTYPE(L).EQ.1) THEN
		   DO 22 J=ISV(L),IEV(L)
		   DO 22 M=1,NCAT(J)
			IF (IIX(II,J,M).EQ.1) THEN
			   PROB(K,II)=THETA(K,J,M)*PROB(K,II)
			END IF
22		   CONTINUE
		ELSE IF (IPARTYPE(L).EQ.3) THEN
		   DO 35 J=ISV(L),IEV(L)
			IF (IVARTYPE(J).EQ.3) THEN
			   DO 36 M=1,NCAT(J)
			     IF (IIX(II,J,M).EQ.1) THEN
			   	PROB(K,II)=THETA(K,J,M)*PROB(K,II)
*			ELSE IF (IIX(II,J,M).EQ.9) THEN
*			   IF (THETA(K,J,M).GT.ALIM
*    :			     PROB(K,II)=(THETA(K,J,M)**THETA(K,J,M))
*     :						*PROB(K,II)
			     END IF
36			   CONTINUE
			END IF
35		   CONTINUE
		END IF
37	      CONTINUE

*	*******HAVEN'T DONE THE MISSING PART FOR THE LOCATION MODEL 
*	*******CHECK CATEGORICAL CONTRIBNTO LIKELIHOOD FOR MISSING***

*	evaluate the continuous variables contribution to the densities

	      DO 23 L=1,NPAR
		IF (IPARTYPE(L).NE.1) THEN
		  DENS(K,L)=0.0
		  I1=0

*	check whether variables are missing in a partition, and create 
*	ISWPCOL(j). The augmented covariance matrix is swept on the 
*	columns specified in ISWPCOL.

		IF (IPARTYPE(L).EQ.2) THEN
		   ISUMISS=0
		   J1=0
		   DO 42 J=ISV(L),IEV(L)
			J1=J1+1
			ISUMISS=ISUMISS + IMISSCODE(II,J)
			ISWPCOL(J1) = IMISSCODE(II,J)
42		   CONTINUE

*	If some variables are missing in the partition, sweep and then 
*	do the regression

		IF ((ISUMISS.GT.0).AND.(ISUMISS.LE.IPC(L))) THEN
     		  CALL SWEEP(AUGCOV,ISWPCOL,ISUMISS,K,L,IPC(L),EN)
		  J3 = 0
		  DO 43 J = ISV(L),IEV(L)
		     J3=J3+1
		     IF (IMISSCODE(II,J).EQ.1) THEN
			ESTX(K,II,J)=EN(K,L,0,J3)
			J2=0
			DO 44 I=ISV(L),IEV(L)
			   J2=J2+1
			   IF (IMISSCODE(II,I).EQ.0) ESTX(K,II,J) = 
     :				ESTX(K,II,J) + EN(K,L,J3,J2)*X(II,I)
44			CONTINUE
		     END IF
43		  CONTINUE
		END IF
		END IF

*	Look at the contribution of the observed values to the likelihood
*	First, the inverse and the determinant of the covariance matrix 
*	for the observed values of observation I is found by using NAG.  

	IF (IPARTYPE(L).EQ.2) THEN
	   I1=0
	   I3=0
	   DO 45 I=ISV(L),IEV(L)
		I3=I3+1
		I2=0
		IF (IMISSCODE(II,I).EQ.0) THEN
		   I1=I1+1
		   I4=0
		   DO 46 J=ISV(L),IEV(L)
			I4=I4+1
			IF(IMISSCODE(II,J).EQ.0) THEN
			   I2=I2+1
			   VAROBS(K,L,I1,I2)=VARIX(K,L,I3,I4)
			END IF
46		   CONTINUE
		END IF
45	   CONTINUE
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
			DO 24 I=ISV(L),IEV(L)
			    J1=0
			    I1=I1+1
			    IF (IMISSCODE(II,I).EQ.0) THEN
				I2=I2+1
				J2=0
			  	DO 47 J=ISV(L),IEV(L)
				   J1=J1+1
				   IF (IMISSCODE(II,J).EQ.0) THEN
					J2=J2+1
	           			DENS(K,L) = DENS(K,L) 
     :    + (X(II,I)-EMU(K,L,I1))*VARIN(K,L,I2,J2)*(X(II,J)-EMU(K,L,J1))
				   END IF
47				CONTINUE
			    END IF
24			  CONTINUE
*	(ii) evaluate the continuous location variables contribution
		  ELSE IF(IPARTYPE(L).EQ.3) THEN
			DO 25 I=ISV(L),IEV(L)
			 IF(IVARTYPE(I).EQ.3) M=IX(II,I)
25			CONTINUE
			DO 26 I=ISV(L),IEV(L)
			 IF (IVARTYPE(I).EQ.4) THEN
			   J1=0
			   I1=I1+1
			   DO 30 J=ISV(L),IEV(L)
			    IF (IVARTYPE(J).EQ.4) THEN
			     J1=J1+1
			     DENS(K,L)=DENS(K,L) + (X(II,I)-EMUL(K,L,I1,M))
     :  *VARIN(K,L,I1,J1)*(X(II,J)-EMUL(K,L,J1,M))
			    END IF
30			   CONTINUE
			 END IF
26			CONTINUE
		  END IF
	IF (ISUMISS.LT.IPC(L)) THEN
		  DENS(K,L) = DEXP(-0.5*DENS(K,L))
		  A=0.5*FLOAT(IPC(L)-ISUMISS)
		  DENS(K,L)=DENS(K,L)/((2.0*PIE)**(A)*DSQRT(ADET(K,L)))
		  APRODENS(K)=APRODENS(K)*DENS(K,L)
	END IF
		END IF
23	     CONTINUE
	     APRODENS(K)=APRODENS(K)*PROB(K,II)*PI(K)
	     SDENS=SDENS+APRODENS(K)
21	   CONTINUE
	   IF (SDENS.NE.0.0) THEN 
		DO 27 K=1,NG
		    Z(II,K)=APRODENS(K)/SDENS
27		CONTINUE
		ALL=ALL+DLOG(SDENS)
	   ELSE
		DO 28 K=1,NG
		   Z(II,K)=0.01
28		CONTINUE
		   WRITE(7,29)II
29		   FORMAT(//' SUM OF DENSITY FUNCTIONS IS ZERO',
     :		   ' FOR OBSERVATION',I4)
	   END IF
20	CONTINUE

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

100 	DO 101 K=1,NG
	 ZSUM(K)=0.0
	   DO 102 II=1,NOBS
	  	ZSUM(K)=ZSUM(K) + Z(II,K)
102	   CONTINUE
	 PI(K) = ZSUM(K)/NOBS
101  	CONTINUE

*	(ii) the conditional probabilities

	DO 103 K=1,NG
	DO 103 L=1,NPAR
		IF(IPARTYPE(L).EQ.1) THEN
		  DO 104 J=ISV(L),IEV(L)
		  DO 104 M=1,NCAT(J)
			  THETA(K,J,M) = 0.0
			  DO 105 II=1,NOBS
			     IF(IMISSCODE(II,J).EQ.1) THEN
				THETA(K,J,M) = THETA(K,J,M) + 
     :				Z(II,K)*THETA2(K,J,M)
			     ELSE IF(IMISSCODE(II,J).EQ.0) THEN
				THETA(K,J,M) = THETA(K,J,M) +
     :				 Z(II,K)*IIX(II,J,M)
			     END IF
105			  CONTINUE
			  THETA(K,J,M)=THETA(K,J,M)/ZSUM(K)
104		  CONTINUE
		ELSE IF(IPARTYPE(L).EQ.2) THEN

*	(iii) the means (EMU(K,L,J))

		  J1=0
		  DO 106 J=ISV(L),IEV(L)
		   J1=J1+1
		   XSUM(K,J1)=0.0
		   DO 107 II=1,NOBS
		      IF(IMISSCODE(II,J).EQ.1) X(II,J)=ESTX(K,II,J)
		      XSUM(K,J1)=XSUM(K,J1) + X(II,J)*Z(II,K)
107		   CONTINUE
		   EMU(K,L,J1)=XSUM(K,J1)/ZSUM(K)
106		  CONTINUE
*	(iv) the location model parameters
*	   (i)the discrete variable

		ELSE IF(IPARTYPE(L).EQ.3) THEN
		  DO 108 J=ISV(L),IEV(L)
		    IF(IVARTYPE(J).EQ.3) THEN
			DO 109 M=1,NCAT(J)
			  THETA(K,J,M)=0.0
			  DO 110 II=1,NOBS
			     IF(IX(II,J).EQ.M) THETA(K,J,M) 
     :					= THETA(K,J,M) + Z(II,K)
110			  CONTINUE
			  THETA(K,J,M)=THETA(K,J,M)/ZSUM(K)
109		        CONTINUE
		    END IF
108		  CONTINUE

*		(ii) the continuous variables the means EMUL(K,L,J,M)
		DO 111 J=ISV(L),IEV(L)
		  IF (IVARTYPE(J).EQ.3) THEN
		     DO 112 M=1,IM(L)
			J1=0
			DO 112 JJ=ISV(L),IEV(L)
			  IF (IVARTYPE(JJ).EQ.4) THEN
			    J1=J1+1
			    XSUM2(K,J1,M)=0.0
		            ZM(K,M)=0.0
			    DO 113 II=1,NOBS
			      IF(IX(II,J).EQ.M) THEN 
				XSUM2(K,J1,M)=XSUM2(K,J1,M) + 
     :					     X(II,JJ)*Z(II,K)
				ZM(K,M)=ZM(K,M)+Z(II,K)
			      END IF
113			     CONTINUE
			  EMUL(K,L,J1,M)=XSUM2(K,J1,M)/ZM(K,M)
			  END IF
112			CONTINUE
		    END IF
111		CONTINUE
		WRITE(7,601)K,((EMUL(K,L,J1,M),M=1,IM(L)),J1=1,IPC(L))
601	FORMAT(/2X,'FOR GROUP',I2,X,'EMUL',10F10.4)
		END IF
103	CONTINUE


*	Calculate updated estimates of the variances
*	(i) the multivariate normal data

	DO 115 K=1,NG
		DO 116 L=1,NPAR
		 IF (IPARTYPE(L).NE.1) THEN
		   DO 117 J=1,IPC(L)
	 	   DO 117 I=1,J
		        VAR(K,L,I,J)=0.0
117		   CONTINUE
		 END IF
116		CONTINUE
		DO 118 II=1,NOBS
		DO 118 L=1,NPAR
		      IF (IPARTYPE(L).EQ.2) THEN
			J1=0
			DO 119 J=ISV(L),IEV(L)
		           IF(IMISSCODE(II,J).EQ.1) X(II,J)=ESTX(K,II,J)
			   I1=0
			   J1=J1+1
			   DO 119 I=ISV(L),J
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
119			CONTINUE
		   END IF
118		   CONTINUE

*	(ii) the continuous location data

		   DO 120 L=1,NPAR
		      IF (IPARTYPE(L).EQ.3) THEN
			DO 121 J=ISV(L),IEV(L)
			 IF (IVARTYPE(J).EQ.3) THEN
			  DO 122 II=1,NOBS
			  DO 122 M=1,IM(L)
			    IF(IX(II,J).EQ.M) THEN
			     J1=0
			     DO 123 JJ=ISV(L),IEV(L)
				IF (IVARTYPE(JJ).EQ.4) THEN
			         I1=0
			         J1=J1+1
			          DO 124 I=ISV(L),J
			           IF (IVARTYPE(I).EQ.4) THEN
			            I1=I1+1
			      VAR(K,L,I1,J1) = VAR(K,L,I1,J1)+(X(II,JJ) 
     :			      -EMUL(K,L,J1,M))*(X(II,I)-EMUL(K,L,I1,M))
     :					*Z(II,K)
			    	   END IF
124				  CONTINUE
				END IF
123			     CONTINUE
			    END IF
122			  CONTINUE
			 END IF
121			CONTINUE
		      END IF
120		   CONTINUE
		DO 126 L=1,NPAR
		 IF (IPARTYPE(L).NE.1) THEN
		   DO 125 J=1,IPC(L)
		   DO 125 I=1,J
			VAR(K,L,I,J)=VAR(K,L,I,J)/ZSUM(K)
			VAR(K,L,J,I)=VAR(K,L,I,J)
125		   CONTINUE
		 END IF
126		CONTINUE
115	CONTINUE

*	Make a copy of the covariance matrix before we use NAG
 
	DO 130 K=1,NG
	DO 130 L=1,NPAR
	    IF(IPARTYPE(L).NE.1) THEN
		DO 131 J=1,IPC(L)
		DO 131 I=1,IPC(L)
		   VARIX(K,L,I,J)=VAR(K,L,I,J)
131	        CONTINUE
	    END IF
130	CONTINUE

*	Make a copy of the conditional probabilities for use with 
*	missing categorical data

	DO 132 K=1,NG
	DO 132 J=1,NVAR
	   IF (IVARTYPE(J).EQ.1) THEN
		DO 133 M=1,NCAT(J)
		   THETA2(K,J,M)=THETA(K,J,M)
133		CONTINUE
	   END IF
132	CONTINUE

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
	DO 507 K=1,NG
	DO 507 L=1,NPAR
	   IF (IPARTYPE(L).EQ.1) THEN
	      DO 508 J=ISV(L),IEV(L)
		WRITE(11,509)(THETA(K,J,M),M=1,NCAT(J))
509		FORMAT(10F10.6)	
508	      CONTINUE
	   ELSE IF(IPARTYPE(L).EQ.2) THEN
	      WRITE(11,510) (EMU(K,L,J),J=1,IPC(L))
510	      FORMAT(10F13.6)
	   ELSE
	      DO 511 J=ISV(L),IEV(L)
		IF (IVARTYPE(J).EQ.3) THEN
		  WRITE(11,509)(THETA(K,J,M),M=1,NCAT(J))
		END IF
511	      CONTINUE
	      DO 512 J=1,IPC(L)
		WRITE(11,510) (EMUL(K,L,J,M),M=1,IM(L))
512	      CONTINUE
	   END IF
507	CONTINUE
	DO 513 K=1,NG
	DO 513 L=1,NPAR
	   IF (IPARTYPE(L).NE.1) THEN
	      DO 514 I=1,IPC(L) 
		 WRITE(11,515)(VARIX(K,L,I,J),J=1,IPC(L))
515		 FORMAT(10F13.6)
514	      CONTINUE
	   END IF
513	CONTINUE

*	(2) the current estimates of the parameters are printed out 
*	Estimates of the proportions in each group and loglikelihood

500	WRITE(7,888),ICNT,CLOGLI(ICNT)
	DO 540 K=1,NG
	   WRITE(7,541) K,PI(K)
541	   FORMAT(/' THE ESTIMATE OF THE MIXING PROPORTION IN GROUP ',I3,
     :	   'IS ', F10.8)
540	CONTINUE

*	estimates of the probabilities for each group

	DO 525 L=1,NPAR
	  WRITE(7,502)
	  WRITE(7,547)L
547	FORMAT(/,' THE CURRENT ESTIMATES FOR PARTITION',I3)
	DO 525 K=1,NG
	    WRITE(7,524)K
524	    FORMAT(/,' Group:',I3,/,' ---------') 
	   IF (IPARTYPE(L).EQ.1) THEN
		DO 528 J=ISV(L),IEV(L)
		  WRITE(7,529)J
529		  FORMAT(' For variable ',I3,X,'the level ',
     :	'probabilities are')
		  WRITE(7,530) (THETA(K,J,M),M=1,NCAT(J))
530		  FORMAT(10F10.6)
528	        CONTINUE
	   ELSE IF(IPARTYPE(L).EQ.2) THEN

*	Estimate of the means for each group.

		WRITE(7,531) (EMU(K,L,J),J=1,IPC(L))
531		FORMAT(' The mean for this partition is ',/,10F13.6)

*	estimates of the location model parameters

	ELSE IF(IPARTYPE(L).EQ.3) THEN
		DO 533 J=ISV(L),IEV(L)
		  IF (IVARTYPE(J).EQ.3) THEN
		   WRITE(7,534)J,(THETA(K,J,M),M=1,NCAT(J))
534		   FORMAT(' For variable',I3,X,'THETA(K,J,M)'
     :					' is',10F10.6)
		  END IF
533		CONTINUE
		J1=0
		DO 535 J=ISV(L),IEV(L)
		 IF(IVARTYPE(J).EQ.4) THEN
		   J1=J1+1
		   WRITE(7,536) (EMUL(K,L,J1,M),M=1,IM(L)) 
536		   FORMAT(' The mean for the continuous location'
     :				' variables is',/,10F13.6)
		 END IF
535		CONTINUE
	END IF
525	CONTINUE

*	Estimates of the variances for each group.

	DO 523 L=1,NPAR
	   IF(IPARTYPE(L).NE.1) THEN
	        WRITE(7,502)
		WRITE(7,522)L
522		FORMAT(/,' THE CURRENT ESTIMATE OF THE COVARIANCE MATRIX FOR',
     :  ' PARTITION ',I3)
	DO 543 K=1,NG
	   WRITE(7,526)K
526	   FORMAT(/,' Group',I3,':',/,' --------')
	   DO 543 I=1,IPC(L)
		   WRITE(7,527)(VAR(K,L,I,J),J=1,IPC(L))
527		   FORMAT(10F15.6)
543	   CONTINUE
	      END IF
523	CONTINUE
	WRITE(7,502)

*	Determine the assignment of the observations to groups

	DO 516 I=1,NOBS
		IMAX = 1
		   DO 517 K=2,NG
		      IF(Z(I,K).GT.Z(I,IMAX)) THEN
		         IMAX=K
		      END IF
517		   CONTINUE
		IGP(I)=IMAX
516	CONTINUE
	WRITE(7,542)
542	FORMAT(/,X,'THE ASSIGNMENT OF OBSERVATIONS TO GROUPS',/)
	WRITE(7,518) (IGP(I),I=1,NOBS)
518	FORMAT(10I3)
	DO 537 K=1,NG
		NUM(K)=0
537	CONTINUE
	DO 538 I=1,NOBS
	DO 538 K=1,NG
		IF(IGP(I).EQ.K) NUM(K)=NUM(K)+1
538	CONTINUE
	WRITE(7,539)(NUM(K),K=1,NG)
539	FORMAT(1X,/,' TOTAL NUMBERS IN EACH GROUP',/,10I5)

*	The estimates of the missing data given the group assignment

	WRITE(7,502)
	WRITE(7,544)
544	FORMAT(/,' THE ESTIMATES OF THE MISSING DATA VALUES GIVEN THE ',
     :  'GROUP ASSIGNMENT',//,
     :	5X,'GROUP',5X,'OBSERVATION',4X,'VARIABLE',3X,'ESTIMATED VALUE')
	DO 545 II=1,NOBS
	DO 545 J=1,NVAR
	DO 545 K=1,NG
	    IF (IMISSCODE(II,J).EQ.1) THEN
	      IF (IGP(II).EQ.K) THEN
		IF (IVARTYPE(J).EQ.2) THEN
		   WRITE(7,546)K,II,J,ESTX(K,II,J)
546		   FORMAT(7X,I2,8X,I5,10X,I2,6X,F16.8)
		END IF
	      END IF
	    END IF
545	CONTINUE
	WRITE(7,502)

*	have a look at the Zij,s, and write out the assigned groups and 
*	the Zij's to GROUPS.OUT

	WRITE(7,521)
521	FORMAT(//' THE ESTIMATES OF THE POSTERIOR PROBALITIES')
		DO 519 I=1,NOBS
		   WRITE(7,520)I,(Z(I,K),K=1,NG)
520		   FORMAT(' OBSERVATION',I4,2X,10F10.6)
		   WRITE(12,532) IGP(I), (Z(I,K),K=1,NG)
532		   FORMAT(X,I2,X,10F9.6)
519		CONTINUE
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
  	    	DO 202 I=1,IDIM
	    	DO 202 J=1,IDIM
		   TEMP(I,J) = VAROBS(IK,IL,I,J)
202	    	CONTINUE
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
		   DO 213 I = 1,IDIM
		   DO 213 J = 1,IDIM
		      VARIN(IK,IL,I,J)=B(I,J)
213		   CONTINUE
		ELSE
		   WRITE(7,207)IFAIL,IK,IL
207		   FORMAT(//,' TO CALCULATE THE INVERSE IFAIL = ',I2, 
     :		   ' FOR GROUP', I2,' PARTITION', I3)
		END IF
	DO 199 J=1,IDIM
	     DO 199 I=1,J
	        VARIN(IK,IL,I,J)=VARIN(IK,IL,J,I)
199	CONTINUE 
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

	DO 315 I=0,IC
	DO 315 J=0,IC
	   EM(IK,IL,I,J)=AUG(IK,IL,I,J)
315	CONTINUE

	IF (ISUM.LT.IC) THEN
	  DO 303 J1=1,IC
	    IF (ISP(J1).EQ.0) THEN
		I=J1

*	Calculation of the diagonal elements

		EN(IK,IL,I,I) = -1/EM(IK,IL,I,I)

*	Remainder of the iTH row and column

		DO 300 J=0,IC
		  IF(J.NE.I) THEN 
			EN(IK,IL,I,J)=EM(IK,IL,I,J)/EM(IK,IL,I,I)
			EN(IK,IL,J,I)=EN(IK,IL,I,J)	
		  END IF
300		CONTINUE

*	The remaining elements of N

		DO 301 J=0,IC
	    	   IF (J.NE.I) THEN
			DO 302 K=0,J
		   	   IF(K.NE.I) THEN
				EN(IK,IL,J,K) = EM(IK,IL,J,K) 
     :  				- EM(IK,IL,J,I) * EN(IK,IL,I,K)
				EN(IK,IL,K,J)=EN(IK,IL,J,K)
		   	   END IF
302			CONTINUE
	    	   END IF
301		CONTINUE
		DO 305 JJ=0,IC
		DO 305 II=0,IC
		    EM(IK,IL,II,JJ)=EN(IK,IL,II,JJ)
305		CONTINUE
	    END IF
303	  CONTINUE
	ELSE IF(ISUM.EQ.IC) THEN
	  DO 304 J=0,IC
	  DO 304 J2=0,J
	      EN(IK,IL,J,J2)=EM(IK,IL,J,J2)
	      EN(IK,IL,J2,J)=EN(IK,IL,J,J2)
304	  CONTINUE
	END IF
	RETURN
	END
