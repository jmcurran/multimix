          	progAM MISSING &
 
!	This program fits a mixture of multivariate distributions 
!	using the EM algorithm to a data set that has data missing at 
!	random. The data file contains both categorical and continuous 
!	variables. The missing data must be coded to the value set in 
!	NOTE (2). 
 
!	If the program does not converge after iter=200 iterations, 
!	the estimates of the parameters will be entered into 
!	EMPARAMEST.OUT. This file can then be used as the parameter 
!	input file for PROGRAM MISSING if desired. 
 
!	The assignment of the observations to groups (IGP(i)) and the 
!	posterior probabilities (Zij's) are entered into GROUPS.OUT. 
!	This file can be used in MINITAB etc. for further analysis. 
 
!	NOTE: 
 
!	(1) THIS PROGRAM REQUIRES VARIABLES IN A PARTITION TO BE STORED 
!	CONTIGUOUSLY. HENCE THE DATA IS READ IN WITH THE VARIABLE ORDER 
!	BEING SPECIFIED BY JP(J). IVARTYPE(J) AND NCAT(J) BOTH REFER TO 
!	THE REARRANGED DATA 
 
!	(2) The code for missing data is -999. imiss=-999. 
 
!	(3) IF SDENS=0, THEN Z(II,K) IS SET TO 0.01 
 
!	(4) THE PROGRAM CURRENTLY HAS A MAXIMUM OF 
!		1500 OBSERVATIONS 		(iob=1500) 
!		6 GROUPS	  		(ik6=6) 
!		15 ATTRIBUTES & 15 PARTITIONS 	(ip15=15) 
!		10 LEVELS OF CATEGORIES 	(im10=10) 
!		200 ITERATIONS FOR CONVERGENCE 	(iter=200) 
!         ******ALTER IF REQUIRED****** 
!	  (REMEMBER TO ALTER PARAMETERS IN DETINV ALSO) 
 
!	The parameter file contains:- 
!	NG - the number of groups 
!	NOBS - the number of observations 
!	NVAR - the number of variables 
!	NPAR - the number of partitions 
!	ISPEC - an indicator variable for a specified grouping of the 
!		observations (1 = observations are not specified into 
!		groups, 2 = observations are specified into groups) 
!	JP(j) - column in the data array in which the jTH variable of 
!	         the file will be stored 
!	IP(l) - number of variables in the lTH partition, l=1,NPAR 
!	IPC(l) - number of continuous variables in partition l, l=1,NPAR 
!	ISV(l) - indicator starting value for the partition, l=1,NPAR 
!	IEV(l) - indicator end value for the partition, l=1,NPAR 
!	IPARTYPE(l) - indicator giving the type of model each partition is 
!		(1 = categorical, 2 = MVN, 3 = location models), l=1,NPAR 
!	IVARTYPE(j) - an indicator variable for type of variable, j=1,NVAR 
!		(1 = categorical variable, 2 = continuous variable, 
!		3 = categorical variable involved in location model, 
!		4 = continuous variable involved in the location model) 
!	NCAT(j) - number of categories of jTH variable. 
!		(For continuous variables, set NCAT(J)=0) 
!	PI(k) - estimated mixing proportions for each group, k=1,NG 
!	THETA(K,J,M) - estimated probability that the jTH categorical 
!			variable is at level M, given that in group k 
!	EMUL(k,l,j,m)- estimated mean vector for each group, each 
!	partition, and each location of the location model variables 
!	EMU(K,L,J) -estimated mean vector for each group and each 
!	partition of continuous variables 
!	VARIX(K,L,I,J)   -estimated covariance matrices for each group 
 
              	impICIT DOUBLE PRECISION(A-H,O-Z) &
              	charCTER*16 infile, datafile &
              	paraETER (PIE=3.141592653589792, iob=1500, ip15=15, ik6 = 6, &
              im10=10, iter=200, imiss = -9) &
              	dimeSION EMU(ik6,ip15,ip15),VAR(ik6,ip15,ip15,ip15), &
              Z(iob,ik6), PI(ik6), DENS(ik6,ip15), VARIN(ik6,ip15,ip15,ip15), &
              ZSUM(ik6), XSUM(ik6,ip15), ADET(ik6,ip15), CLOGLI(iter), &
              VARIX(ik6,ip15,ip15,ip15), APRODENS(ik6), IP(ip15), IPC(ip15), &
              ISV(ip15), IEV(ip15), IPARTYPE(ip15), EMUL(ik6,ip15,ip15,im10), &
              THETA(ik6,ip15,0:im10), PROB(ik6,iob),IVARTYPE(ip15),NCAT(ip15), &
              IM(ip15), IGP(iob), NUM(ik6), XSUM2(ik6,ip15,ip15),ZM(ik6,ip15), &
              IGRP(iob), JP(ip15), IX(iob,ip15), X(iob,ip15), &
              IMISSCODE(iob,ip15), IIX(iob,ip15,im10), ISWPCOL(ip15), &
              	AUGCOV(ik6,ip15,0:ip15,0:ip15),EN(ik6,ip15,0:ip15,0:ip15), &
              VAROBS(ik6,ip15,ip15,ip15), ESTX(ik6,iob,ip15), &
              THETA2(ik6,ip15,im10) &
              	prin *, ' MIXTURE ESTIMATION BY EM' &
              	prin *, '-------------------------------' &
              	prin *, 'Data file: ' 
    	read(*,10) datafile &
              10	foMAT (A16) &
              	prin *, 'Parameter file: ' 
    	read(*,10) infile 
    	open(7, FILE='GENERAL.OUT', STATUS = 'NEW') 
    	open(8, FILE=datafile, STATUS = 'OLD') 
    	open(9, FILE=infile, STATUS = 'OLD') &
              	open12, FILE='GROUPS.OUT',STATUS = 'NEW') 
    	read(9,*) NG, NOBS, NVAR, NPAR,ISPEC 
    	read(9,*) (JP(J),J = 1,NVAR) 
    	read(9,*) (IP(L),L = 1,NPAR) 
    	read(9,*) (IPC(L), L = 1,NPAR) 
    	read(9,*) (ISV(L),L = 1,NPAR) 
    	read(9,*) (IEV(L),L = 1,NPAR) 
    	read(9,*) (IPARTYPE(L),L = 1,NPAR) 
    	read(9,*) (IVARTYPE(J),J = 1,NVAR) 
    	read(9,*) (NCAT(J),J = 1,NVAR) &
 
!	read in estimates of the parameters if grouping not specified 
 
              	if(ipec == 1) then &
              	  reD(9,*) (PI(K),K = 1,NG) 
    	  do12 K = 1,NG 
    	  do12 L = 1,NPAR &
              		ifipartype(l) == 1) then 
    			do 13 J = ISV(L),IEV(L) &
              EAD(9,*) (THETA(K,J,M),M = 1,NCAT(J)) &
              13			ONTINUE &
              		els if(ipartype(l) == 2) then &
              			reD(9,*)(EMU(K,L,J),J = 1,IPC(L)) &
              		els if(ipartype(l) == 3) then 
    			do14 J = ISV(L),IEV(L) &
              f(ivartype(j) == 3) then &
              				rAD(9,*)(THETA(K,J,M),M = 1,NCAT(J)) &
              				i(L) = NCAT(J) 
    END IF &
              14			ONTINUE &
              			j10 
    			do15 J = ISV(L),IEV(L) &
              f(ivartype(j) == 4) then 
    J1 = J1+1 
    READ(9,*) &
              			    (EMUL(K,L,J1,M),M = 1,IM(L)) &
              ND IF &
              15			ONTINUE 
    		endIF &
              12	coTINUE &
 
              	do 1 K = 1,NG &
              	do 1 L = 1,NPAR &
              		if(partype(l) /= 1) &
              		READ(9,*)((VARIX(K,L,I,J),J=1,IPC(L)),I = 1,IPC(L)) &
              16  	ONTINUE &
 
    	else 
              	  reD(9,*) (IGRP(I),I = 1,NOBS) 
    	  do18 I = 1,NOBS 
    	  do18 K = 1,NG 
    Z(I,K) = 0.0 &
              18	ONTINUE 
    	  do31 I = 1,NOBS &
              K = IGRP(I) &
              (I,IK) = 1.0 &
              31	ONTINUE 
    	  do33 L = 1,NPAR 
    if(ipartype(l) == 3) then &
                		do2 J = ISV(L),IEV(L) &
                		  i (ivartype(j) == 3) im(l) = ncat(j) &
                32		cNTINUE &
                ND IF &
                33	ONTINUE &
                	endF &
 
!	read in the observations 
 
                	do 1 I = 1,NOBS &
                		rea(8,* ) (X(I,JP(J)),J = 1,NVAR) &
                17	coTINUE &
 
!	Check for the presence of missing values, and create 
!	IMISSCODE(I,J), where IMISSCODE(I,J)=1 if observation i variable j is 
!	missing, and IMISSCODE(I,J)=0 otherwise 
 
                	do 3 I = 1,NOBS &
                	do 3 J = 1,NVAR &
                		iva=NINT(X(I,J)) &
                		ifival == imiss) then &
                			imSSCODE(I,J) = 1 &
                		els &
                			imSSCODE(I,J) = 0 
      		endIF &
                34	coTINUE 
 
!	Check the categorical variable to see if 0 is a category. 
!	Add 1 to the variables if it is so that Xij=0 is taken as 
!	level 1, Xij=1 is level 2 etc. 
 
      	do 8 J = 1,NVAR &
                f((ivartype(j) == 1).or.(ivartype(j) == 3)) then &
                		imi=NINT(X(1,J)) &
                		ifimin /= 0) then &
                O 81 I = 2,NOBS &
                			inXT = NINT(X(I,J)) &
                			ifinext == 0) imin = inext 
    end do 
    		endIF &
              		ifimin == 0) then 
    do i = 1,nobs 
      			if(imisscode(i,j) == 0) then 
      IX(I,J) = NINT(X(I,J)) + 1 &
                			elE 
      IX(I,J) = NINT(X(I,J)) &
                			en IF 
    end do &
              		els 
    do i = 1,nobs &
                			ixI,J) = NINT(X(I,J)) 
    end do 
    		endIF &
              	   eD IF &
              80	coTINUE &
 
!	For the categorical variables, create the matrix IIX(I,J,M) where 
!	x(i,j,m)=1 if obs. i, variable j is at level m, 
!	x(i,j,m)=0 if not at level m. 
 
              	do 8 I = 1,NOBS &
              	do 8 J = 1,NVAR &
              f((ivartype(j) == 1).or.(ivartype(j) == 3)) then &
              		do5 M = 1,NCAT(J) 
    if (ix(i,j) == m) then &
                			ii(I,J,M) = 1 
    else if (ix(i,j) == imiss) then &
                			ii(I,J,M) = 9 
    ELSE &
                			ii(I,J,M) = 0 
    END IF &
              85		cNTINUE 
    END IF &
              84	coTINUE &
 
!	print the input values 
 
              	writ(7,50)NG,NOBS,NVAR,NPAR &
              50	foMAT(/,' NO OF GROUPS IS ',I3,/,' NO OF OBSERVATIONS IS ', &
              I5,/,' NO OF VARIABLES ',I3,/,' NO OF PARTITIONS', I3) &
              	writ(7,51)(IP(L),L = 1,NPAR) &
              51	foMAT(/,' THE NO OF VARIABLES IN EACH PARTITION IS',/,10I3) &
 
!	Send to the M step if ISPEC = 2 (groups specified) 
              	if(ipec == 2) go to 100 &
!	otherwise, print the parameter estimates 
              	writ(7,52) &
              52	foMAT(/,' MIXING PROPORTIONS') &
              	writ(7,53)(PI(K),K = 1,NG) &
              53	foMAT(/,10F6.3) &
              	do 5 K = 1,NG &
              	do 5 L = 1,NPAR &
              		if(partype(l) == 1) then &
              		  wITE(7,54) K 
    54		FORMAT(/,2X,' THETA(K,J,M) FOR GROUP ',I2) &
              		  d 64 J = ISV(L),IEV(L) 
    WRITE(7,55)(THETA(K,J,M),M = 1,NCAT(J)) 
    55		    FORMAT(10F8.4) 
    end do &
              		els if(ipartype(l) == 2) then &
              		  wITE(7,57) K,L 
    57		FORMAT(/,' FOR GROUP ',I2, X, ' AND PARTITION',I2,X, &
              		  ' THE MEAN IS') &
              		  wITE(7,58)(EMU(K,L,J),J = 1,IPC(L)) 
    58		FORMAT(10F12.6) &
              		els &
              		  wITE(7,57) K,L &
              		  wITE(7,58)((EMUL(K,L,J,M),M=1,IM(L)),J = 1,IPC(L)) &
              		  d 62 J = ISV(L),IEV(L) &
              f (ivartype(j) == 3) then 
    WRITE(7,54) 
    WRITE(7,55)(THETA(K,J,M),M = 1,NCAT(J)) &
              ND IF 
    end do &
              	   eD IF &
              56	coTINUE &
              	do 5 K = 1,NG &
              	do 5 L = 1,NPAR 
    if (ipartype(l) /= 1) then &
                		wriE(7,60) L,K &
                60		fRMAT(/,2X,'VARIANCE FOR PARTITION',I2,X,' AND GROUP',I2) &
                		do3 I = 1,IPC(L) &
                RITE(7,61)(VARIX(K,L,I,J),J = 1,IPC(L)) 
      61		 FORMAT(10F13.6) &
                63		cNTINUE 
    END IF &
              59	coTINUE &
 
!	Print out the types of models for each partition 
 
              	do 7 L = 1,NPAR &
              	   i (ipartype(l) == 1) then &
              		wriE(7,75) L,IP(L),IPARTYPE(L) &
              75		fRMAT(' PARTITION',I3,' HAS',I2,' VARIABLES',/,' ITYPE', &
              ' IS',I2,' HENCE A CATEGORICAL MODEL FOR THIS PARTITION') &
              	   ese if (ipartype(l) == 2) then &
              		wriE(7,76) L,IP(L),IPARTYPE(L) &
              76		fRMAT(' PARTITION',I3,' HAS',I2,' VARIABLES',/,' ITYPE', &
              ' IS',I2,' HENCE A MVN MODEL FOR THIS PARTITION') &
              	   eSE &
              		wriE(7,77) L,IP(L),IPARTYPE(L) &
              77		fRMAT(' PARTITION',I3,' HAS',I2,'VARIABLES',/,' ITYPE', &
              ' IS',I2,' HENCE A LOCATION MODEL FOR THIS PARTITION') &
              	endF &
              74	coTINUE &
 
!	E step of the EM algorithm. 
!	Estimate the complete data sufficient statistics given the 
!	data, & current values of means,variances and mixing proportions. 
 
              	icnt1 
 
!	Form the augmented covariance matrix 
 
    99999do k = 1,ng 
      	do 4 L = 1,NPAR &
                f (ipartype(l) == 2) then &
                		augOV(K,L,0,0) = -1. &
                		do1 J = 1,IPC(L) 
      AUGCOV(K,L,J,0) = EMU(K,L,J) 
      AUGCOV(K,L,0,J) = AUGCOV(K,L,J,0) 
      do j1 = 1,j 
        AUGCOV(K,L,J,J1) = VARIX(K,L,J,J1) 
        AUGCOV(K,L,J1,J) = AUGCOV(K,L,J,J1) &
                  41		cNTINUE &
                  ND IF &
                  40	coTINUE 
 
 
        	all=.0 &
                  	alim1E-7 
        	do 2 II = 1,NOBS &
                  	   sENS = 0.0 &
                  	   d 21 K = 1,NG 
        PROB(K,II) = 1.0 
        APRODENS(K) = 1.0 
 
!	evaluate the discrete variables contribution to the densities 
 
        do l = 1,npar &
                    		ifipartype(l) == 1) then &
                    O 22 J = ISV(L),IEV(L) &
                    O 22 M = 1,NCAT(J) 
          			if(iix(ii,j,m) == 1) then 
          PROB(K,II) = THETA(K,J,M)*PROB(K,II) &
                    			en IF 
        end do &
                  		els if (ipartype(l) == 3) then &
                  O 35 J = ISV(L),IEV(L) 
        			if(ivartype(j) == 3) then 
        do m = 1,ncat(j) 
          if (iix(ii,j,m) == 1) then 
            	PROB(K,II) = THETA(K,J,M)*PROB(K,II) 
!			ELSE IF (IIX(II,J,M).EQ.9) THEN 
!			   IF (THETA(K,J,M).GT.ALIM 
!    :			     PROB(K,II)=(THETA(K,J,M)**THETA(K,J,M)) 
!     :						*PROB(K,II) 
          END IF 
        end do &
                  			en IF 
      end do 
      		endIF 
    end do 
 
!	*******HAVEN'T DONE THE MISSING PART FOR THE LOCATION MODEL 
!	*******CHECK CATEGORICAL CONTRIBNTO LIKELIHOOD FOR MISSING*** 
 
!	evaluate the continuous variables contribution to the densities 
 
    do l = 1,npar &
                		ifipartype(l) /= 1) then &
                		  dNS(K,L) = 0.0 &
                		  i=0 &
 
!	check whether variables are missing in a partition, and create 
!	ISWPCOL(j). The augmented covariance matrix is swept on the 
!	columns specified in ISWPCOL. 
 
                		ifipartype(l) == 2) then &
                SUMISS = 0 &
                1 = 0 &
                O 42 J = ISV(L),IEV(L) &
                			j1J1+1 &
                			isMISS = ISUMISS + IMISSCODE(II,J) &
                			isPCOL(J1) = IMISSCODE(II,J) 
    end do &
 
!	If some variables are missing in the partition, sweep and then 
!	do the regression 
 
              		if(isumiss > 0).and.(isumiss <= ipc(l))) then 
    	  CALL SWEEP(AUGCOV,ISWPCOL,ISUMISS,K,L,IPC(L),EN) &
              		  j = 0 &
              		  d 43 J = ISV(L),IEV(L) 
    J3 = J3+1 
    if (imisscode(ii,j) == 1) then &
                			esX(K,II,J) = EN(K,L,0,J3) &
                			j20 
      			do44 I = ISV(L),IEV(L) 
      J2 = J2+1 
      if (imisscode(ii,i) == 0) estx(k,ii,j) = &
                				ESTX(K,II,J) + EN(K,L,J3,J2)*X(II,I) &
                44			ONTINUE 
    END IF 
    end do 
    		endIF 
    		endIF &
 
!	Look at the contribution of the observed values to the likelihood 
!	First, the inverse and the determinant of the covariance matrix 
!	for the observed values of observation I is found by using NAG. 
 
              	if (partype(l) == 2) then &
              	   i=0 &
              	   i=0 &
              	   d 45 I = ISV(L),IEV(L) &
              		i3=3+1 
    		i2= &
              		ifimisscode(ii,i) == 0) then &
              1 = I1+1 &
              4 = 0 &
              O 46 J = ISV(L),IEV(L) &
              			i4I4+1 &
              			ifimisscode(ii,j) == 0) then 
    I2 = I2+1 
    VAROBS(K,L,I1,I2) = VARIX(K,L,I3,I4) &
              			en IF 
    end do 
    		endIF 
    end do &
              	ndimIPC(L)-ISUMISS &
              	if (sumiss == 0) then &
              	   cLL DETINV(VAROBS,K,L,IPC(L),ADET,VARIN) 
    	elseif (isumiss < ipc(l)) then &
                	   cLL DETINV(VAROBS,K,L,NDIM,ADET,VARIN) &
                	endF &
                	endF &
 
!	(i) evaluate the MVN contribution 
                		  i(ipartype(l) == 2) then &
                			i10 &
                			i20 
      			do24 I = ISV(L),IEV(L) 
      J1 = 0 
      I1 = I1+1 
      if (imisscode(ii,i) == 0) then &
                  				i=I2+1 &
                  				j=0 
        do j = isv(l),iev(l) 
          J1 = J1+1 
          if (imisscode(ii,j) == 0) then &
                      2 = J2+1 
            			DENS(K,L) = DENS(K,L) &
                      + (X(II,I)-EMU(K,L,I1))*VARIN(K,L,I2,J2)*(X(II,J)-EMU(K,L,J1)) 
          END IF 
        end do 
      END IF 
    end do &
!	(ii) evaluate the continuous location variables contribution 
              		  ese if(ipartype(l) == 3) then 
    			do25 I = ISV(L),IEV(L) &
              			 i(ivartype(i) == 3) m = ix(ii,i) &
              25			ONTINUE 
    			do26 I = ISV(L),IEV(L) &
              			 i (ivartype(i) == 4) then 
    J1 = 0 
    I1 = I1+1 
    do j = isv(l),iev(l) 
      if (ivartype(j) == 4) then 
        J1 = J1+1 
        DENS(K,L) = DENS(K,L) + (X(II,I)-EMUL(K,L,I1,M)) &
                  *VARIN(K,L,I1,J1)*(X(II,J)-EMUL(K,L,J1,M)) 
      END IF 
    end do &
              			 eD IF &
              26			ONTINUE &
              		  eD IF &
              	if (sumiss < ipc(l)) then &
              		  dNS(K,L) = DEXP(-0.5*DENS(K,L)) &
              		  a0.5*FLOAT(IPC(L)-ISUMISS) &
              		  dNS(K,L) = DENS(K,L)/((2.0*PIE)**(A)*DSQRT(ADET(K,L))) &
              		  aRODENS(K) = APRODENS(K)*DENS(K,L) &
              	endF 
    		endIF 
    end do 
    APRODENS(K) = APRODENS(K)*PROB(K,II)*PI(K) 
    SDENS = SDENS+APRODENS(K) 
    end do &
              	   i (sdens /= 0.0) then &
              		do7 K = 1,NG 
    Z(II,K) = APRODENS(K)/SDENS &
              27		cNTINUE &
              		allALL+DLOG(SDENS) &
              	   eSE &
              		do8 K = 1,NG &
              (II,K) = 0.01 &
              28		cNTINUE &
              RITE(7,29)II 
    29		 FORMAT(//' SUM OF DENSITY FUNCTIONS IS ZERO', &
              		   ' FOR OBSERVATION',I4) &
              	   eD IF &
              20	coTINUE &
 
!	Check on convergence - look at the likelihood function 
!	If the absolute value of the difference in 2 likelihoods is 
!	less than a tolerance value the estimates are written out. 
!	A check is made on the number of iterations (200 max). 
!	statement 100 - the M step 
!	statement 500 - print out the current estimates (algorithm 
!	has converged) 
!	statement 501 - algorithm hasn't converged, current estimates 
!	are printed out. 
 
              	clogI(ICNT) = ALL &
              	writ(7,888),ICNT,CLOGLI(ICNT) &
              888	fRMAT(1X,'FOR LOOP',I5,'LOGLIKELIHOOD IS',F16.8) &
              	if (cnt <= 10) go to 100 &
              	c=1.-10 &
              	tol=BS(CLOGLI(ICNT) - CLOGLI(ICNT-10)) &
              	if(tl <= c) go to 500 &
              	if(int >= iter) go to 501 &
 
!	M step of the EM algorithm 
!	Calculate an updated estimate of the parameters 
!	 means and covariances. 
 
!	(i) the mixing proportions, 
 
              100 	O 101 K = 1,NG &
              	 zsu(K) = 0.0 &
              	   d 102 II = 1,NOBS &
              	  	zUM(K) = ZSUM(K) + Z(II,K) 
    end do &
              	 pi() = ZSUM(K)/NOBS 
    end do 
 
!	(ii) the conditional probabilities 
 
    	do 13 K = 1,NG 
    	do 13 L = 1,NPAR &
              		if(partype(l) == 1) then &
              		  d 104 J = ISV(L),IEV(L) &
              		  d 104 M = 1,NCAT(J) &
              HETA(K,J,M) = 0.0 &
              O 105 II = 1,NOBS 
    if(imisscode(ii,j) == 1) then &
                				tETA(K,J,M) = THETA(K,J,M) + &
                				Z(II,K)*THETA2(K,J,M) 
    else if(imisscode(ii,j) == 0) then &
                				tETA(K,J,M) = THETA(K,J,M) + &
                				 Z(II,K)*IIX(II,J,M) 
    END IF 
    end do &
              HETA(K,J,M) = THETA(K,J,M)/ZSUM(K) 
    end do &
              		els if(ipartype(l) == 2) then &
 
!	(iii) the means (EMU(K,L,J)) 
 
              		  j=0 &
              		  d 106 J = ISV(L),IEV(L) &
              1 = J1+1 &
              SUM(K,J1) = 0.0 &
              O 107 II = 1,NOBS 
    if(imisscode(ii,j) == 1) x(ii,j) = estx(k,ii,j) 
    XSUM(K,J1) = XSUM(K,J1) + X(II,J)*Z(II,K) 
    end do &
              MU(K,L,J1) = XSUM(K,J1)/ZSUM(K) 
    end do &
!	(iv) the location model parameters 
!	   (i)the discrete variable 
 
              		els if(ipartype(l) == 3) then &
              		  d 108 J = ISV(L),IEV(L) 
    if(ivartype(j) == 3) then 
      			do109 M = 1,NCAT(J) &
                HETA(K,J,M) = 0.0 &
                O 110 II = 1,NOBS 
      if(ix(ii,j) == m) theta(k,j,m) &
                					= THETA(K,J,M) + Z(II,K) 
    end do &
              HETA(K,J,M) = THETA(K,J,M)/ZSUM(K) 
    end do 
    END IF 
    end do &
 
!		(ii) the continuous variables the means EMUL(K,L,J,M) 
              		do11 J = ISV(L),IEV(L) &
              		  i (ivartype(j) == 3) then 
    do m = 1,im(l) &
                			j10 
      			do112 JJ = ISV(L),IEV(L) &
                f (ivartype(jj) == 4) then 
      J1 = J1+1 
      XSUM2(K,J1,M) = 0.0 
      ZM(K,M) = 0.0 
      do ii = 1,nobs 
        if(ix(ii,j) == m) then &
                    				xUM2(K,J1,M) = XSUM2(K,J1,M) + &
                    					     X(II,JJ)*Z(II,K) &
                    				z(K,M) = ZM(K,M)+Z(II,K) 
        END IF 
      end do &
                MUL(K,L,J1,M) = XSUM2(K,J1,M)/ZM(K,M) &
                ND IF 
    end do 
    END IF &
              111		ONTINUE &
              		wriE(7,601)K,((EMUL(K,L,J1,M),M=1,IM(L)),J1 = 1,IPC(L)) &
              601	fRMAT(/2X,'FOR GROUP',I2,X,'EMUL',10F10.4) 
    		endIF &
              103	cNTINUE &
 
 
!	Calculate updated estimates of the variances 
!	(i) the multivariate normal data 
 
              	do 15 K = 1,NG &
              		do16 L = 1,NPAR 
    		 if(ipartype(l) /= 1) then &
              O 117 J = 1,IPC(L) 
    do i = 1,j 
      VAR(K,L,I,J) = 0.0 
    end do &
              		 en IF &
              116		ONTINUE &
              		do18 II = 1,NOBS &
              		do18 L = 1,NPAR 
    if (ipartype(l) == 2) then &
                			j10 
      			do119 J = ISV(L),IEV(L) 
      if(imisscode(ii,j) == 1) x(ii,j) = estx(k,ii,j) 
      I1 = 0 
      J1 = J1+1 
      do i = isv(l),j 
        	      if(imisscode(ii,i) == 1)x(ii,i) = estx(k,ii,i) 
        I1 = I1+1 
        if ((imisscode(ii,i) == 1) .and. &
                  				(imisscode(ii,j) == 1)) then &
                  				vR(K,L,I1,J1) = VAR(K,L,I1,J1) + &
                  ((X(II,J) - EMU(K,L,J1))*(X(II,I)-EMU(K,L,I1)) &
                  + EN(K,L,I1,J1))*Z(II,K) 
      ELSE 
        VAR(K,L,I1,J1) = VAR(K,L,I1,J1) + (X(II,J) &
                  			        -EMU(K,L,J1))*(X(II,I)-EMU(K,L,I1))*Z(II,K) 
      END IF 
    end do &
              ND IF 
    end do &
 
!	(ii) the continuous location data 
 
              O 120 L = 1,NPAR 
    if (ipartype(l) == 3) then 
      			do121 J = ISV(L),IEV(L) &
                			 i (ivartype(j) == 3) then &
                O 122 II = 1,NOBS &
                O 122 M = 1,IM(L) 
      if(ix(ii,j) == m) then 
        J1 = 0 
        do jj = isv(l),iev(l) &
                    				i (ivartype(jj) == 4) then 
          I1 = 0 
          J1 = J1+1 
          do i = isv(l),j 
            if (ivartype(i) == 4) then 
              I1 = I1+1 
              VAR(K,L,I1,J1) = VAR(K,L,I1,J1)+(X(II,JJ) &
                        			      -EMUL(K,L,J1,M))*(X(II,I)-EMUL(K,L,I1,M)) &
                        					*Z(II,K) 
              	   END IF 
              	  end do &
                        				eD IF 
            end do 
          END IF 
        end do &
                  			 eD IF 
      end do 
    END IF 
    end do &
              		do26 L = 1,NPAR 
    		 if(ipartype(l) /= 1) then &
              O 125 J = 1,IPC(L) &
              O 125 I = 1,J &
              			va(K,L,I,J) = VAR(K,L,I,J)/ZSUM(K) &
              			va(K,L,J,I) = VAR(K,L,I,J) 
    end do &
              		 en IF &
              126		ONTINUE &
              115	cNTINUE &
 
!	Make a copy of the covariance matrix before we use NAG 
 
              	do 10 K = 1,NG &
              	do 10 L = 1,NPAR &
              f(ipartype(l) /= 1) then &
              		do31 J = 1,IPC(L) &
              		do31 I = 1,IPC(L) &
              ARIX(K,L,I,J) = VAR(K,L,I,J) 
    end do &
              ND IF &
              130	cNTINUE &
 
!	Make a copy of the conditional probabilities for use with 
!	missing categorical data 
 
              	do 12 K = 1,NG &
              	do 12 J = 1,NVAR &
              	   i (ivartype(j) == 1) then &
              		do33 M = 1,NCAT(J) &
              HETA2(K,J,M) = THETA(K,J,M) &
              133		ONTINUE &
              	   eD IF &
              132	cNTINUE &
 
              	icntICNT+1 &
 
!	send back to the E step 
              	go t 99999 &
 
!	Write out the current estimates of the parameters. 
 
!	(1) If the algorithm has not converged:- 
 
              501	wITE(7,502) &
              502	fRMAT(/'----------------------------------------------------') &
              	writ(7,503) iter &
              503	fRMAT(/' THE EM ALGORITHM HAS NOT CONVERGED AFTER ',I3, &
              /,' ITERATIONS BUT THE CURRENT ESTIMATES OF THE PARAMETERS ', &
              	/,' WILL BE PRINTED OUT.') &
              	writ(7,502) &
 
!	The parameters are to be written to 'EMPARAMEST.DAT' to be 
!	used as input for the PROGRAM MULTIMIX. ISPEC is set to 1. 
 
              	open11, FILE='EMPARAMEST.OUT',STATUS = 'NEW') &
              	ispe=1 &
              	writ(11,504) NG, NOBS, NVAR, NPAR, ISPEC &
              504	fRMAT(X,5I6) &
              	writ(11,505) (JP(J),J = 1,NVAR) &
              	writ(11,505) (IP(L),L = 1,NPAR) &
              	writ(11,505) (IPC(L),L = 1,NPAR) &
              	writ(11,505) (ISV(L),L = 1,NPAR) &
              	writ(11,505) (IEV(L),L = 1,NPAR) &
              	writ(11,505) (IPARTYPE(L),L = 1,NPAR) &
              	writ(11,505) (IVARTYPE(J),J = 1,NVAR) &
              	writ(11,505) (NCAT(J),J = 1,NVAR) &
              505	fRMAT(10I4) &
              	writ(11,506) (PI(K),K = 1,NG) &
              506	fRMAT(10F10.6) 
    	do 57 K = 1,NG 
    	do 57 L = 1,NPAR &
              	   i (ipartype(l) == 1) then 
    do j = isv(l),iev(l) &
                		wriE(11,509)(THETA(K,J,M),M = 1,NCAT(J)) &
                509		ORMAT(10F10.6) 
    end do &
              	   ese if(ipartype(l) == 2) then 
    WRITE(11,510) (EMU(K,L,J),J = 1,IPC(L)) 
    510	    FORMAT(10F13.6) &
              	   eSE 
    do j = isv(l),iev(l) &
                		ifivartype(j) == 3) then &
                		  wITE(11,509)(THETA(K,J,M),M = 1,NCAT(J)) 
      		endIF 
    end do 
    do j = 1,ipc(l) &
                		wriE(11,510) (EMUL(K,L,J,M),M = 1,IM(L)) 
    end do &
              	   eD IF &
              507	cNTINUE &
              	do 53 K = 1,NG &
              	do 53 L = 1,NPAR &
              	   i (ipartype(l) /= 1) then 
    do i = 1,ipc(l) &
                		 wrTE(11,515)(VARIX(K,L,I,J),J = 1,IPC(L)) 
      515		FORMAT(10F13.6) 
    end do &
              	   eD IF &
              513	cNTINUE &
 
!	(2) the current estimates of the parameters are printed out 
!	Estimates of the proportions in each group and loglikelihood 
 
              500	wITE(7,888),ICNT,CLOGLI(ICNT) &
              	do 50 K = 1,NG &
              	   wITE(7,541) K,PI(K) 
    541	 FORMAT(/' THE ESTIMATE OF THE MIXING PROPORTION IN GROUP ',I3, &
              	   'IS ', F10.8) &
              540	cNTINUE &
 
!	estimates of the probabilities for each group 
 
              	do 55 L = 1,NPAR &
              	  wrTE(7,502) &
              	  wrTE(7,547)L &
              547	fRMAT(/,' THE CURRENT ESTIMATES FOR PARTITION',I3) &
              	do 55 K = 1,NG &
              RITE(7,524)K 
    524	  FORMAT(/,' Group:',I3,/,' ---------') &
              	   i (ipartype(l) == 1) then &
              		do28 J = ISV(L),IEV(L) &
              		  wITE(7,529)J 
    529		 FORMAT(' For variable ',I3,X,'the level ', &
              	'probabilities are') &
              		  wITE(7,530) (THETA(K,J,M),M = 1,NCAT(J)) 
    530		 FORMAT(10F10.6) 
    end do &
              	   ese if(ipartype(l) == 2) then &
 
!	Estimate of the means for each group. 
 
              		wriE(7,531) (EMU(K,L,J),J = 1,IPC(L)) &
              531		ORMAT(' The mean for this partition is ',/,10F13.6) 
 
!	estimates of the location model parameters 
 
    	elseif(ipartype(l) == 3) then &
                		do33 J = ISV(L),IEV(L) &
                		  i (ivartype(j) == 3) then &
                RITE(7,534)J,(THETA(K,J,M),M = 1,NCAT(J)) 
      534		  FORMAT(' For variable',I3,X,'THETA(K,J,M)' &
                					' is',10F10.6) &
                		  eD IF &
                533		ONTINUE 
      		j1= &
                		do35 J = ISV(L),IEV(L) &
                		 ifivartype(j) == 4) then &
                1 = J1+1 &
                RITE(7,536) (EMUL(K,L,J1,M),M = 1,IM(L)) 
      536		  FORMAT(' The mean for the continuous location' &
                				' variables is',/,10F13.6) &
                		 en IF &
                535		ONTINUE &
                	endF &
                525	cNTINUE &
 
!	Estimates of the variances for each group. 
 
                	do 53 L = 1,NPAR &
                	   i(ipartype(l) /= 1) then 
      WRITE(7,502) &
                		wriE(7,522)L &
                522		ORMAT(/,' THE CURRENT ESTIMATE OF THE COVARIANCE MATRIX FOR', &
                ' PARTITION ',I3) &
                	do 53 K = 1,NG &
                	   wITE(7,526)K 
      526	 FORMAT(/,' Group',I3,':',/,' --------') &
                	   d 543 I = 1,IPC(L) &
                RITE(7,527)(VAR(K,L,I,J),J = 1,IPC(L)) 
      527		  FORMAT(10F15.6) 
    end do 
    END IF &
              523	cNTINUE &
              	writ(7,502) &
 
!	Determine the assignment of the observations to groups 
 
              	do 56 I = 1,NOBS &
              		ima = 1 &
              O 517 K = 2,NG 
    if(z(i,k) > z(i,imax)) then 
      IMAX = K 
    END IF 
    end do &
              		igpI) = IMAX &
              516	cNTINUE &
              	writ(7,542) &
              542	fRMAT(/,X,'THE ASSIGNMENT OF OBSERVATIONS TO GROUPS',/) &
              	writ(7,518) (IGP(I),I = 1,NOBS) &
              518	fRMAT(10I3) &
              	do 57 K = 1,NG &
              		numK) = 0 &
              537	cNTINUE &
              	do 58 I = 1,NOBS &
              	do 58 K = 1,NG &
              		if(gp(i) == k) num(k) = num(k)+1 &
              538	cNTINUE &
              	writ(7,539)(NUM(K),K = 1,NG) &
              539	fRMAT(1X,/,' TOTAL NUMBERS IN EACH GROUP',/,10I5) &
 
!	The estimates of the missing data given the group assignment 
 
              	writ(7,502) &
              	writ(7,544) &
              544	fRMAT(/,' THE ESTIMATES OF THE MISSING DATA VALUES GIVEN THE ', &
              'GROUP ASSIGNMENT',//, &
              	5X,'GROUP',5X,'OBSERVATION',4X,'VARIABLE',3X,'ESTIMATED VALUE') &
              	do 55 II = 1,NOBS &
              	do 55 J = 1,NVAR &
              	do 55 K = 1,NG &
              f (imisscode(ii,j) == 1) then 
    if (igp(ii) == k) then &
                		ifivartype(j) == 2) then &
                RITE(7,546)K,II,J,ESTX(K,II,J) 
      546		  FORMAT(7X,I2,8X,I5,10X,I2,6X,F16.8) 
      		endIF 
    END IF &
              ND IF &
              545	cNTINUE &
              	writ(7,502) &
 
!	have a look at the Zij,s, and write out the assigned groups and 
!	the Zij's to GROUPS.OUT 
 
              	writ(7,521) &
              521	fRMAT(//' THE ESTIMATES OF THE POSTERIOR PROBALITIES') &
              		do19 I = 1,NOBS &
              RITE(7,520)I,(Z(I,K),K = 1,NG) 
    520		  FORMAT(' OBSERVATION',I4,2X,10F10.6) &
              RITE(12,532) IGP(I), (Z(I,K),K = 1,NG) 
    532		  FORMAT(X,I2,X,10F9.6) &
              519		ONTINUE &
    	end 
 
!	This subroutine calculates the inverse and the determinant of 
!	the covariance matrix using NAG routines F03ABF and F01ABF. 
 
              	subrUTINE DETINV(VAROBS,IK,IL,IDIM,ADET,VARIN) &
              	implCIT DOUBLE PRECISION(A-H,O-Z) &
              	paraETER (ik6=6,ip15=15,IIP = ip15+1) &
              	dimeSION VAROBS(ik6,ip15,ip15,ip15),VARIN(ik6,ip15,ip15,ip15), &
              ADET(ik6,ip15) &
              	inteER IA,IFAIL,IB &
              	real8 TEMP(IIP,IIP),DET,WKSPCE(IIP),ZI(IIP),B(IIP,IIP) &
              	ia=iP &
              	ib=iP &
              	ifai = 0 
    N = IDIM 
    	DO 202 I = 1,IDIM 
    do j = 1,idim &
                EMP(I,J) = VAROBS(IK,IL,I,J) 
      	end do &
                		cal F03ABF(TEMP,IA,N,DET,WKSPCE,IFAIL) &
                		conINUE &
                		ififail == 0) then &
                DET(IK,IL) = DET &
                		els &
                RITE(7,204)IFAIL,IK,IL 
      204		  FORMAT(//,' TO CALCULATE THE DETERMINANT IFAIL = ',I2, &
                		   ' FOR GROUP',I3,' PARTITION',I3) 
      		endIF &
                		ifaL = 0 &
                		cal F01ABF(TEMP,IA,N,B,IB,ZI,IFAIL) &
                		conINUE &
                		if(fail == 0) then &
                O 213 I = 1,IDIM &
                O 213 J = 1,IDIM 
      VARIN(IK,IL,I,J) = B(I,J) 
    end do &
              		els &
              RITE(7,207)IFAIL,IK,IL 
    207		  FORMAT(//,' TO CALCULATE THE INVERSE IFAIL = ',I2, &
              		   ' FOR GROUP', I2,' PARTITION', I3) 
    		endIF &
              	do 19 J = 1,IDIM 
    do i = 1,j 
      VARIN(IK,IL,I,J) = VARIN(IK,IL,J,I) &
                199	cNTINUE &
                	retuN &
      	end 
 
!	This subroutine sweeps the augmented covariance matrix on the 
!	observed values. Columns are specified in ISWPCOL 
!	EN = SWP [i1,i2,...,it] EM where i1,i2,..,it are the columns 
!	specified in ISWPCOL 
 
                	subrUTINE SWEEP(AUG,ISP,ISUM,IK,IL,IC,EN) &
                	implCIT DOUBLE PRECISION(A-H,O-Z) &
                	paraETER(ik6=6,ip15 = 15) &
                	dimeSION EM(ik6,ip15,0:ip15,0:ip15),EN(ik6,ip15,0:ip15,0:ip15), &
                		ISP(ip15),AUG(ik6,ip15,0:ip15,0:ip15) &
 
!	First we make a copy of the augmented covariance matrix 
 
                	do 35 I = 0,IC &
                	do 35 J = 0,IC &
                	   e(IK,IL,I,J) = AUG(IK,IL,I,J) &
                315	cNTINUE &
 
                	if (sum < ic) then 
      	  do303 J1 = 1,IC &
                f (isp(j1) == 0) then &
                		i=j &
 
!	Calculation of the diagonal elements 
 
                		en(K,IL,I,I) = -1/EM(IK,IL,I,I) &
 
!	Remainder of the iTH row and column 
 
                		do00 J = 0,IC &
                		  i(j /= i) then &
                			enIK,IL,I,J) = EM(IK,IL,I,J)/EM(IK,IL,I,I) &
                			enIK,IL,J,I) = EN(IK,IL,I,J) &
                		  eD IF &
                300		ONTINUE &
 
!	The remaining elements of N 
 
                		do01 J = 0,IC 
      if (j /= i) then 
        			do302 K = 0,J 
        if(k /= i) then &
                    				e(IK,IL,J,K) = EM(IK,IL,J,K) &
                    				- EM(IK,IL,J,I) * EN(IK,IL,I,K) &
                    				e(IK,IL,K,J) = EN(IK,IL,J,K) 
        END IF 
      end do 
    END IF &
              301		ONTINUE &
              		do05 JJ = 0,IC &
              		do05 II = 0,IC 
    EM(IK,IL,II,JJ) = EN(IK,IL,II,JJ) &
              305		ONTINUE &
              ND IF 
    end do 
    	elseif(isum == ic) then 
      	  do304 J = 0,IC 
      	  do304 J2 = 0,J 
      EN(IK,IL,J,J2) = EM(IK,IL,J,J2) 
      EN(IK,IL,J2,J) = EN(IK,IL,J,J2) 
    end do &
              	endF &
              	retuN 
    	end 
