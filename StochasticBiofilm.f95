! A fortran95 program for G95
!Mathematical Model of Biofilms
!Part of Doctoral Thesis for Curtin University and University of Aberdeen, completed in February 2021.
! Created by Dr Paulina A Dzianach under supervision and guidance of Dr Francisco Perez-Reche, Prof. Gary A Dykes, Prof. Kenneth Forbes and Prof. Norval Strachan.

!Parameters Legend
!*inistate* : Initial inoculum. Can take values of "single cell" or "inoculation". "single cell" starts the simulation with one cell attached at the centre of the surface. "inoulation" starts the simulation with a given number of live cells spread randomly across the surface, the number specified by the value of inicellcount
!*L, H*: length and width of the system in terms of number of sites
!*maxstep*: maximum number of time steps in the simulation
!*max_movie_step*: maximum number of frames in generated movies of simulated biofilm formation
!*max_rls*: number of realisations for a given parameter set
!*inicellcount* : initial number of cells inoculated onto the surface
!*t_random* : value 1 if the length of time steps is stochastic and 0 if it's deterministic
!*tmax*: maximum time of the simulation (h)
!*gmax*:maximum growth rate of cells (1/h)
!*Kc*:Monod growth saturation coefficient for c uptake and for growth rate equation (mmol/L here converted to [10^-7 mmol/site])
!*Ko*:Monod growth rate sat.coeff for oxygen uptake (mmol/L here converted to [10^-7 mmol/site])
!*umin_c*:nutrient uptake starvation threshold (mmol/h here converted to [10^-7 mmol/(site*h)])
!*c*: carbon source concentration in the bulk medium (mmol/L here converted to [10^-7 mmol/site])
!*o*:oxygen concentration in bulk medium (mmol/L here converted to [10^-7 mmol/site])
!*Al*:lysis A coeff (L/(mmol*h) here converted to [10^7 site/(mmol*h)])
!*Betty*:lysis B coeff (L/(mmol*h) here converted to [10^7 site/(mmol*h)])	
!*dbs_cj*:starvation death rate (1/h)
!*diff_noeps*:diffusion of nutrient/oxygen molecules through the biofilm in absence of ecm (m^2/h)
!*diff_eps*:diffusion of nutrient/oxygen molecules through the biofilm with ecm (m^2/h)
!*dx*: length of a single patch (m)
!*UmaxO*:maximum uptake rate of O molecules in a single patch (mmol/h here converted to [10^-7 mmol/(site*h)])
!*UmaxC*:!maximum uptake rate of C molecules in a single patch (mmol/h here converted to [10^-7 mmol/(site*h)])
!*EPS_max*:maximum EPS extrection rate (1/h)

program main

!!!!CHOOSE MODEL MODE---------ECM uptake/no ECM uptake------------
    character (len=*), parameter :: ECM_uptake="YES" !Valid values: "YES" or "NO"
!!!!
!----PARAMETERS-------------- 
    character (len=*), parameter ::inistate="single cell"   
    integer, parameter :: L=20, H=50,maxstep=10000, max_movie_step=0, max_rls=100, inicellcount=1, t_random=1
    integer, parameter :: tmax=12, NEupt=4		!time of the simulation (h), number of times ECM can be utilized 		 
    real, parameter  :: gmax=0.8, Kc=6.5E-2,Ko=3E-5,umin_c=0.01 		        
    real, parameter :: c=6.5E-2,o=0.265E-2,Al=2,Betty=1000, dbs_cj=0.3			 	        	 	        
    real, parameter :: diff_noeps=9E-10 	
    real*4, parameter :: diff_eps=9E-11,  dx=5E-6 	!diff_eps=9E-10 - faster diffusion, diff_eps=9E-11 - regular diffusion	
    real, parameter :: UmaxO=4.5E-3, UmaxC=4.5E-2					      	
    real, parameter :: EPS_max=0.0 	
    real,parameter :: bonus=0.5!bonus to growth rate from ECM uptake (rE)
!-----VARIABLES-------------------------
    integer :: i, ycor,xcor,counter, dir,dir1,tstep,rls,k,fno, j	!Iterators/coordinates
    integer :: Ncj_att,Ncj						!Cell numbers storage
    integer :: shove, shovu, shovnext,shovnow, signx, signy     	!Shoving algorithm path variables
    integer :: xmin, xmax, endpath, xdis, ydis, movex           	!Shoving algorithm path variables c.d.
    integer :: space, test                                      	!Boolean variables
    integer :: d,r							!Matrix coordinates
    real    :: dentry, rentry,EPSnow, EPSnext				!Matrix element values, uptake of compounds
    integer :: summit,biosurf,inv_surf,newlast                          !summit stores the top of the biofilm
    integer :: last, next, m, pore, top, no_pores, b_area					!Surface enlargement - biosurf determination queue algorithm 
    real    :: SurfEn, biofilm_density
    real    :: Q, Qnext, Q1, Q2, Q3, Q4, dt				!Events: Q is the total rate of all events and Qnext the random number from 0 to Q which determines what exactly will happen
    real    :: eventcount						!Used in all event routines to find a cell which was chosen to act by Qnext and proceed.					
    real    :: oxy, conc 						!Diff : diffusion coef., dependant on biofilm ECM concentration, oxy-local oxygen concentration and conc - local c concentration
    integer, dimension(L*H)   :: checked, queue                         !Surfac enlargement - biosurf determination queue algorithm 
    real, dimension(L*H)      :: gcj, EPS, dbs,dbl,uptk_c,uptk_o	!local rates of events: growth, ECM excretion, deactivation, lysis 
    integer, dimension(2,4) :: neighbor
    integer, dimension(4) :: neighbor_node
    real,dimension(L*H):: concentration, Alpha,Alpha_o, concentration_o,no_EPS
    integer, dimension(L*H) :: cj_attached,Nei, cj_deactivated !used to track a moving particle, cj_attached: attached particles
    integer, dimension(4,L*H) :: finalnnodes
    integer,dimension(L*H,L*H):: concentration_matrix, Mat
    real,dimension(L*H,L*H) :: Matrix, M_inv
    real,dimension(max_movie_step,L*H)::biofilm_frame
    real, dimension(maxstep) ::t,Diff
    real, dimension(max_rls,maxstep) ::rls_t
    integer, dimension(max_rls,maxstep):: rls_bcellcount, rls_eps
    integer, dimension(max_rls,maxstep):: rls_ev_loc
    real, dimension(max_rls,maxstep):: i_uptk_c, i_uptk_o,i_conc, i_cono,mean_uptk_c, mean_uptk_o 
    real, dimension(max_rls,maxstep):: i_growth,i_lysis, i_starvation, mean_c, mean_o
    real, dimension(max_rls,maxstep):: mean_growth,mean_lysis,mean_starvation, waiting_t
    integer, dimension(max_rls):: finalstep
    integer, dimension(max_rls, L*H)   :: rls_invaded                   !All rls - This saves all locations which were ever invaded by bacteria
    real, dimension(max_rls,L*H) ::rls_biopic     !All rls - biofilm at the end of the rls, and compound concentrations
    character(len=3), dimension(max_rls,maxstep):: rls_events !mean_uptk is the mean uptake between all live cells. i uptk is the uptake at location where an event took place

! subroutine EPS generation is outdated and would need to be updated before use. I've introduced a safety net - if it happens that by mistake EPS rate has been changed to non zero, the code will spit out an error
    if (EPS_max>0)then
    	write(*,*) "Please change value of the EPS_max parameter to 0. The subroutine for active EPS generation has not been updated."
        stop
    end if
    call srand(43) !Starts the sequence of random numbers from the 43rd element.
            rls_t=0
            rls_bcellcount=0
            finalstep=0
            rls_events=''
            rls_ev_loc=0
            rls_eps=0
            biofilm_frame=0
	    rls_invaded=0
	    rls_biopic=0
	    i_conc=0
	    i_cono=0
            i_uptk_c=0
            i_uptk_o=0
            mean_uptk_c=0
            mean_uptk_o=0
            mean_c=0
            mean_o=0
            mean_growth=0
            mean_lysis=0
            mean_starvation=0
            do rls=1,max_rls
                write(*,*) rls !displays the current realisation number
                t=0 !starts realisation at time 0
                tstep=1 !starts realisation at step 1
                rls_bcellcount(rls,:)=0 !counts biofilm live cells
                call initialization
                do while (tstep<maxstep) !We will walk for maxstep number of steps
		    !biofilm_frame stores an image of the biofilm at a point in time. Different values are given to different biofilm elements to produce the image
                    if(tstep<=max_movie_step .and. rls==max_rls)then
                    	do i=1,L*H
			    if(cj_attached(i)>0)then
				biofilm_frame(tstep,i)=1
			    else if(no_EPS(i)>0)then
				biofilm_frame(tstep,i)=0.5
			    else if(cj_deactivated(i)>0)then
				biofilm_frame(tstep,i)=0.25
			    end if
                        end do
                    end if
                    
                  !  call motile_bioprox_count !counting the number of motile cells close to the biofilm at every tstep
                    if(ECM_uptake=="NO")then
                    	call Qnext_choosing !Finding Qnext based on the weight of each event
		    else if(ECM_uptake=="YES")then
		    	call Qnext_choosing_ECM_uptake_YES
                    else
			write(*,*) "Please choose a valid mode for the program. Type 'YES' for ECM uptake enabled or 'NO for ECM uptake disabled."
                    end if
                    
		    mean_uptk_c(rls,tstep)=0
		    mean_uptk_o(rls,tstep)=0
                    do i=1,L*H
			if(cj_attached(i)>0)then
				mean_uptk_c(rls,tstep)=mean_uptk_c(rls,tstep)+uptk_c(i) !sum the uptake from all live cells prior to the event taking place at tstep
                        	mean_uptk_o(rls,tstep)=mean_uptk_o(rls,tstep)+uptk_o(i) 
                           	mean_c(rls,tstep)=mean_c(rls,tstep)+concentration(i)
                                mean_o(rls,tstep)=mean_o(rls,tstep)+concentration_o(i)
			end if
                    end do
                    if(Ncj_att>0)then
	            	mean_uptk_c(rls,tstep)=mean_uptk_c(rls,tstep)/Ncj_att !mean uptake of C among live cells
                    	mean_uptk_o(rls,tstep)=mean_uptk_o(rls,tstep)/Ncj_att !mean uptake of O among live cells
                    	mean_c(rls,tstep)=mean_c(rls,tstep)/Ncj_att !mean concentration c among live cells
                    	mean_o(rls,tstep)=mean_o(rls,tstep)/Ncj_att ! mean concentration o among live cells
                   	mean_growth(rls,tstep)=sum(gcj)/Ncj_att !mean growth rate before an event takes place
                   	mean_lysis(rls,tstep)=sum(dbl)/Ncj_att !mean lysis rate before an event takes place
                   	mean_starvation(rls,tstep)=sum(dbs)/Ncj_att !mean starvation rate before an event takes place
                   end if
                   
		  
!--------ALL THE EVENTS------------------------------------------------------------------------------------
                    if (Qnext/=0)then
                        if(Qnext<Q1)then !C.jejuni diffusion in air
                             call cj_growth
			     call finding_summit ! finds the maximum height of the biofilm
                             rls_events(rls,tstep)="Gro"
                             write(*,*) "Growth"
                        else if(Qnext<Q2 .and. Qnext>Q1)then !C.jejuni diffusion in biofilm
                             call EPS_generation
                             call finding_summit
                             rls_events(rls,tstep)="Eps" 
                             write(*,*) "EPS"
                        else if (Qnext>Q2 .and. Qnext<Q3)then !C.jejuni surface attachment
                             call death_by_starvation
                             call finding_summit
                             rls_events(rls,tstep)="Dbs"
                             write(*,*) "Starvation"
                        else if (Qnext>Q3 .and. Qnext<Q4)then
                            call autolysis
                            call finding_summit
                            rls_events(rls,tstep)="Dbl"
                            write(*,*) "Lysis"

                        end if
                        i_uptk_c(rls,tstep)=uptk_c(rls_ev_loc(rls,tstep)) !uptake rates at event location
		  	i_uptk_o(rls,tstep)=uptk_o(rls_ev_loc(rls,tstep))
		        i_conc(rls,tstep)=concentration(rls_ev_loc(rls,tstep)) !chemical concentrations at event location
		        i_cono(rls,tstep)=concentration_o(rls_ev_loc(rls,tstep))
                        i_growth(rls,tstep)=gcj(rls_ev_loc(rls,tstep)) !event rates at event location
                        i_lysis(rls,tstep)=dbl(rls_ev_loc(rls,tstep))
                        i_starvation(rls,tstep)=dbs(rls_ev_loc(rls,tstep))

                    end if
!-------END OF EVENTS----------------------------------------------------
!-----Time incrementation--------------


                  if(tstep/=maxstep)then
                        call realtime_increment
                        call SE_Calc
                        call nutrient_distribution
                        
                     
			rls_bcellcount(rls,tstep)=Ncj_att
                        rls_eps(rls,tstep)=sum(no_EPS)
                        finalstep(rls)=tstep
                        if(t(tstep)>tmax .or. Ncj==0)then
                            tstep=maxstep
                            write(*,*) "reason 1"
                        end if
                        if(tstep>100 .and. Ncj_att==0)then
                            tstep=maxstep
                            write(*,*) "reason 2"
                        end if
			write(*,*) t(tstep), Ncj_att
                        tstep=tstep+1
                    end if

                end do
		
		do i=1,L*H
		    if(cj_attached(i)>0)then
			rls_biopic(rls,i)=1
		    else if(no_EPS(i)>0)then
			rls_biopic(rls,i)=0.5
		    else if(cj_deactivated(i)>1)then
			rls_biopic(rls,i)=0.25
		    end if
		end do
                
		
                rls_t(rls,:)=t(:)

            end do


!__________________________________________________________________________________________________
!                                   DATA STORAGE(RESULTS COLLECTION)
!__________________________________________________________________________________________________
     write(*,*) "collecting data..."
        fno=1
   !     open(fno,file='time_surf_dev_0_3_c_0_3_o_0_3.dat')
    !        do rls=1,max_rls
     !           do tstep=1,finalstep(rls)
      !              write(fno,*) rls_t(rls,tstep), rls
       !         end do
        !    end do
       ! close(fno)
        !fno=fno+1
        !if(fno==6)then
        !	fno=fno+1
       ! end if
     !   open(3,file='GNUbiopic_c2_o10_tmax24_L50H20_10rls_Kcj8_Ko5_dbs3_lys3_Da12_De4.dat')
      !      do i=1,L*H
	!	if(mod(i,L)/=0)then
         !       	write(3,*) mod(i,L),int(i/L)+1,cj_attached(i)+0.5*no_EPS(i) !x,y,z
	!	else
	!		write(3,*) L,int(i/L),cj_attached(i)+0.5*no_EPS(i)+0.25*cj_deactivated(i) ! x,y,z
	!	end if
         !   end do
       ! close(3)



        !end if_
       open(fno,file='fbiopic_v3_1_s1_c_6_5_Aero_tmax_12_rls_100_diff_E11_rE_0_5_betty_1000_Kg_6_5.csv')
 	    do rls=1,max_rls
		do i=1,L*H           
                    write(fno,*) rls_biopic(rls,i), rls
		end do
            end do
        close(fno)
        fno=fno+1
        if(fno==6)then
       		fno=fno+1
        end if



        open(fno,file='data_v3_1_s1_c_6_5_Aero_tmax_12_rls_100_diff_E11_rE_0_5_betty_1000_Kg_6_5.csv')
            do rls=1,max_rls
                do tstep=1,finalstep(rls)
                    write(fno,*) rls_bcellcount(rls,tstep),',', rls,',', rls_eps(rls,tstep),',',&
                                &rls_events(rls,tstep),',', rls_ev_loc(rls,tstep),',',rls_t(rls,tstep),',',&
                                &mean_uptk_c(rls,tstep),',', mean_uptk_o(rls,tstep),',',&
                                &i_uptk_c(rls,tstep),',',i_uptk_o(rls,tstep),',',&
                                &i_growth(rls,tstep),',',i_lysis(rls,tstep),',',&
                                &i_starvation(rls,tstep),',',i_conc(rls,tstep),',',&
                                &i_cono(rls,tstep),',', mean_c(rls,tstep),',',&
                                &mean_o(rls,tstep),',', Diff(tstep),',',&
                                &mean_growth(rls,tstep),',',mean_lysis(rls,tstep),',',&
                                &mean_starvation(rls,tstep),',',waiting_t(rls,tstep)
                end do
            end do
        close(fno)
        fno=fno+1
        if(fno==6)then
        	fno=fno+1
        end if









!__________________________________________________________________________________________________________________________
!                                        ALL THE SUBROUTINES DEFINED BELOW
!             THE PROGRAM DOESN'T DO ANYTING BELOW THIS POINT, IT IS JUST A STORAGE OF DEFINITIONS USED EARLIER
!__________________________________________________________________________________________________________________________
contains
!-----------------SYSTEM INITIALIZATION------------------------------------------------------------
    subroutine initialization
        i=0
        do ycor=1,H !we fix y (the row number) and move on with the column number (x) to make sure the nodes are numbered correctly
            do xcor=1,L
                i=i+1 !this defines our node number as we move on with coordinates
    !General rules to find neighbors of a node with coordinates (xcor,ycor)

                neighbor(1,1)=xcor
                neighbor(2,1)=ycor+1

                neighbor(1,2)=xcor+1
                neighbor(2,2)=ycor

                neighbor(1,3)=xcor
                neighbor(2,3)=ycor-1

                neighbor(1,4)=xcor-1
                neighbor(2,4)=ycor

!---------Applying Boundary conditions (4 closed walls)--------------------------------------------------
                if (xcor==1)then
                    neighbor(1,4)=0
                    neighbor(2,4)=1
                end if
                if(xcor==L)then
                    neighbor(1,2)=0
                    neighbor(2,2)=1
                end if
                if (ycor==1)then
                    neighbor(2,3)=1 !(closed)
                    neighbor(1,3)=0
                end if
                if (ycor==H) then
                    neighbor(2,1)=1
                    neighbor(1,1)=0
                end if
!for top and bottom boundaries, the finalnnodes of the 'neighbours' outside the boundary will return 0.
                neighbor_node(:)=(neighbor(2,:)-1)*L+neighbor(1,:)
                finalnnodes(:, i)=neighbor_node

            end do
        end do
!-------End of applying boundary conditions--------------------------------------------
!------Number of neighbours of each node (size of neighbourhood)-----------------------
        Nei=0
        do i=1,L*H
            do k=1,4
                if(finalnnodes(k,i)>0)then
                    Nei(i)=Nei(i)+1
                end if
            end do
        end do
!------Neighbourhood Matrix M for concentration calculation
        Mat=0 !most of the elements will be 0
        do i=1,L*H
            Mat(i,i)=(-1)*Nei(i)
            do k=1,4
                if(finalnnodes(k,i)/=0)then
                    Mat(i,finalnnodes(k,i))=1
                end if
            end do
        end do
!------End of constructing concentration matrix M


!--------------------------------------------------------------------------------------
        cj_attached=0
        if(inistate=="inoculation")then
            Ncj_att=0
            if(inicellcount>L)then
                do i=1,L
			cj_attached(i)=1
		end do
                Ncj_att=L
		do while(Ncj_att<inicellcount)
			i=0
			do while(i==0)
                    		i=nint(L*rand())
               		end do
			if(i<=L)then
				i=i+L
			end if
               		if(cj_attached(i)<1)then
                   		 cj_attached(i)=cj_attached(i)+1
                	end if
                	Ncj_att=sum(cj_attached)
		end do

            else
            	do while(Ncj_att<inicellcount)
                	 i=0
               		do while(i==0)
                    		i=nint(L*rand())
                	end do
                	if(cj_attached(i)<1)then
                   		 cj_attached(i)=cj_attached(i)+1
                	end if
                	Ncj_att=sum(cj_attached)
           	end do
            end if
	end if
        if(inistate=="single cell")then
        	Ncj_att=1
                i=nint(0.5*L)
                cj_attached(i)=1
        end if

        Ncj=Ncj_att
        no_EPS=0 !initially, no EPS
        cj_deactivated=0 !Initially, all cells are active
        concentration=c !initially, concentration of nutrient is c everywhere
        concentration_o=o !initially, concentration of oxygen is o everywhere
        uptk_c=UmaxC*c/(c+Kc) !total initial uptake of c
        uptk_o=UmaxO*o/(o+Ko) !total initial uptake of o
        summit=2*L

    end subroutine initialization



 !------------------C.JEJUNI GROWTH------------------------------------------------------------------------------------------
    subroutine cj_growth       
	i=0
        eventcount=0
        do while(eventcount<Qnext .and. i<L*H )
            i=i+1
            eventcount=eventcount+gcj(i) 
        end do
        rls_ev_loc(rls,tstep)=i
        Ncj_att=Ncj_att+1
        Ncj=Ncj+1
        !To save time, let's first check if the neighbours of i are empty. If so, applying the shoving algorithm will be unnecessary
        test=0
        counter=0
        dir=int(rand()*4) !dir is an integer from 0 to 3
        do while(test==0 .and. counter<4) ! we will check 4 neighbours, unless test=1, that is unless we find a vacancy earlier
            dir=dir+1 !now dir can be anything from 1 to 4 initially,and it is incremented by 1 each time afterwards
            counter=counter+1 !keep track of the number of neigbours checked for a vacancy
            if(dir>4)then
                dir=dir-4  !make sure we stay within 1-4 range, as we have 4 neigbhours only (at most)
            end if
            if(finalnnodes(dir,i)>0)then ! we must make sure that the neighbour exists (in case i is a boundary site)
                j=finalnnodes(dir,i)
                if(cj_attached(j)+no_EPS(j)==0)then !if the neighbour is empty, we will put the daughter cell there.
                   test=1
                   cj_attached(j)=1
		   i=j
                end if
            end if
        end do
        !FIRST WE FIND THE RANDOM SITE ON THE SURFACE OF THE BIOFILM FOR THE DIRECTED PATH END POINT
        !Finding the perimeter of the microcolony surrounding the cell
        if(test==0)then !in case the neihbours of i don't have a vacancy, we must apply the shoving algorithm
            xmin=i
            xmax=i
!find the x-range of the microcolony to which i belongs (xmax and xmin)
            if(mod(i,L)/=0)then
                do while(cj_attached(xmin)+no_EPS(xmin)>0 .and. xmin>int(i/L)*L+1) !it will stop when a boundary of biofilm or the system is reached
                    xmin=xmin-1
                end do
                do while(cj_attached(xmax)+no_EPS(xmax)>0 .and. xmax<int(i/L)*L+L)
                    xmax=xmax+1
                end do
            else
                do while(cj_attached(xmin)+no_EPS(xmin)>0 .and. xmin>i-L+1)
                    xmin=xmin-1
                end do
           !xmax is just i if i is on the RHS boundary, no brainer.
            end if
!normalize xmax and xmin to numbers between 1 and L
            k=xmax-xmin !k is at most L-1
            
	    xmin=mod(xmin,L) !starting from the surface. xmin is never divisible by L by default. This is at most L-1
            xmax=xmin+k
        !find random point between xmin and xmax
            endpath=xmin+int(rand()*(k+1))
        !go up until vacancy is found
            do while(cj_attached(endpath)+no_EPS(endpath)>0)
                endpath=endpath+L
                if(endpath>L*H)then
                	write(*,*) "Out of the system, out of luck"
                        stop
                end if
            end do
            if(cj_attached(endpath)>0)then
		write(*,*) "endpath wrong"
		stop
	    end if
        !find the distance in xcor and ycor from i to endpath
            signx=1
            signy=1
            xdis=0
            ydis=int(endpath/L)-int(i/L)
	    if(mod(endpath,L)==0 .and. mod(i,L)/=0)then
                ydis=ydis-1
	    else if(mod(endpath,L)/=0 .and. mod(i,L)==0)then
		ydis=ydis+1
	    end if
            if(mod(endpath,L)/=0 .and. mod(i,L)/=0)then
                xdis=mod(endpath,L)-mod(i,L) !this can be negative, or 0
            else if(mod(endpath,L)==0 .and. mod(i,L)/=0)then
                xdis=L-mod(i,L)
            else if(mod(endpath,L)/=0 .and.mod(i,L)==0)then
                xdis=mod(endpath,L)-L
            end if

            if(xdis<0)then
                xdis=abs(xdis)
                signx=-1
            end if

            if(ydis<0)then
                ydis=abs(ydis)
                signy=-1
            end if
            shovnow=1 !first we shove a cell into the site in the path
            EPSnow=0
            do while(xdis>0 .or. ydis>0)
                movex=0
                if(rand()>0.5 .or. ydis==0)then
                    movex=1
                    if(xdis==0)then
                        movex=0
                    end if
                end if
                if(movex==1)then !we shove 1 step in the x direction
                   i=i+signx
                   if(i>L*H)then
                   	write(*,*) "got out of borders x"
                        stop

                   end if
           !see what is already there
                    if(cj_attached(i)+no_EPS(i)==0)then !If we stumble upon a vacancy on the way, just fill the vacancy and end shoving
                        if(shovnow==1)then
                        	cj_attached(i)=1
                        else 
                                no_EPS(i)=EPSnow
                        end if
                        xdis=1
                        ydis=0
                        endpath=i
                    else if(no_EPS(i)>0 .and.shovnow==1)then!_________C--->E___________
                        shovnext=0!next step we'll be shoving EPS
                        EPSnext=no_EPS(i) !this is how much EPS we're going to move next round
                        no_EPS(i)=0
                        cj_attached(i)=1
                    else if(no_EPS(i)>0 .and. shovnow==0)then!_______E--->E____________
                        shovnext=0!next step we'll be shoving EPS
                        if(no_EPS(i)+EPSnow<=1)then !surplus has been distributed, can stop shoving
                                no_EPS(i)=no_EPS(i)+EPSnow
                        	xdis=1
                        	ydis=0
                        	endpath=i
                        else if (no_EPS(i)+EPSnow>1)then
                                no_EPS(i)=1	
                        	EPSnext=(no_EPS(i)+EPSnow)-1 !this is how much EPS we're going to move next round
                                
                        end if
                    else if(cj_attached(i)>0 .and. shovnow==1)then!_____C--->C___________
                        shovnext=1
                        cj_attached(i)=1
                    else if(cj_attached(i)>0 .and. shovnow==0)then!_____E--->C__________
                        shovnext=1
                    	cj_attached(i)=0
                        no_EPS(i)=EPSnow

                    end if
                    xdis=xdis-1
               
              else if(movex==0)then !we shove 1 step in the y direction
                    i=i+signy*L
                   if(i>L*H)then
                   	write(*,*) "got out of borders y"
                        stop

                   end if
                   if(cj_attached(i)+no_EPS(i)==0)then !If we stumble upon a vacancy on the way, just fill the vacancy and end shoving
                        if(shovnow==1)then
                        	cj_attached(i)=1
                        else 
                                no_EPS(i)=EPSnow
                        end if
                        xdis=0
                        ydis=1
                        endpath=i
                    else if(no_EPS(i)>0 .and.shovnow==1)then!_________C--->E___________
                        shovnext=0!next step we'll be shoving EPS
                        EPSnext=no_EPS(i) !this is how much EPS we're going to move next round
                        no_EPS(i)=0
                        cj_attached(i)=1
                    else if(no_EPS(i)>0 .and. shovnow==0)then!_______E--->E____________
                        shovnext=0!next step we'll be shoving EPS
                        if(no_EPS(i)+EPSnow<=1)then !surplus has been distributed, can stop shoving
                                no_EPS(i)=no_EPS(i)+EPSnow
                        	xdis=0
                        	ydis=1
                        	endpath=i
                        else if (no_EPS(i)+EPSnow>1)then
                                no_EPS(i)=1	
                        	EPSnext=(no_EPS(i)+EPSnow)-1 !this is how much EPS we're going to move next round
                                
                        end if
                    else if(cj_attached(i)>0 .and. shovnow==1)then!_____C--->C___________
                        shovnext=1
                        cj_attached(i)=1

                    else if(cj_attached(i)>0 .and. shovnow==0)then!_____E--->C__________
                        shovnext=1
                    	cj_attached(i)=0
                        no_EPS(i)=EPSnow
                    end if
                    if(cj_attached(i)+no_EPS(i)>1)then
                        write(*,*) "Unlucky. Tight squeeze", "attached",cj_attached(i),"EPS", no_EPS(i), i
                        stop
                    end if               
                    ydis=ydis-1
                end if
                shovnow=shovnext  
                EPSnow=EPSnext 
            end do
        end if
	rls_invaded(rls,i)=1 !we save the invaded location
	if(Ncj_att/=sum(cj_attached))then
        	write(*,*) "boohoo", Ncj_att, sum(cj_attached),shovnext, EPSnext, movex, ydis
                stop
        end if


    end subroutine cj_growth
 
   subroutine EPS_generation
        i=0
        eventcount=Q1
        do while(eventcount<Qnext .and. i<L*H)
            i=i+1
            eventcount=eventcount+EPS(i)
        end do
        rls_ev_loc(rls,tstep)=i

        !To save time, let's first check if the neighbours of i are empty. If so, applying the shoving algorithm will be unnecessary
        test=0
        counter=0
        dir=int(rand()*4) !dir is an integer from 0 to 3
        do while(test==0 .and. counter<4) ! we will check 4 neighbours, unless test=1, that is unless we find a vacancy earlier
            dir=dir+1 !now dir can be anything from 1 to 4 initially,and it is incremented by 1 each time afterwards
            counter=counter+1 !keep track of the number of neigbours checked for a vacancy
            if(dir>4)then
                dir=dir-4  !make sure we stay within 1-4 range, as we have 4 neigbhours only (at most)
            end if
            if(finalnnodes(dir,i)>0)then ! we must make sure that the neighbour exists (in case i is a boundary site)
                j=finalnnodes(dir,i)
                if(cj_attached(j)+no_EPS(j)==0)then !if the neighbour is empty, we will put the EPS particle there.
                   test=1
                   no_EPS(j)=1
		   i=j
                end if
            end if
        end do
        !FIRST WE FIND THE RANDOM SITE ON THE SURFACE OF THE BIOFILM FOR THE DIRECTED PATH END POINT
        !Finding the perimeter of the microcolony surrounding the cell
        if(test==0)then !in case the neihbours of i don't have a vacancy, we must apply the shoving algorithm
            xmin=i
            xmax=i
!find the x-range of the microcolony to which i belongs (xmax and xmin)
            if(mod(i,L)/=0)then
                do while(cj_attached(xmin)+no_EPS(xmin)>0 .and. xmin>int(i/L)*L+1) !it will stop when a boundary of biofilm or the system is reached
                    xmin=xmin-1
                end do
                do while(cj_attached(xmax)+no_EPS(xmax)>0 .and. xmax<int(i/L)*L+L)
                    xmax=xmax+1
                end do
            else
                do while(cj_attached(xmin)+no_EPS(xmin)>0 .and. xmin>i-L+1)
                    xmin=xmin-1
                end do
           !xmax is just i if i is on the RHS boundary, no brainer.
            end if
!normalize xmax and xmin to numbers between 1 and L
            k=xmax-xmin !k is at most L-1
            
	    xmin=mod(xmin,L) !starting from the surface. xmin is never divisible by L by default. This is at most L-1
            xmax=xmin+k
        !find random point between xmin and xmax
            endpath=xmin+int(rand()*(k+1))
        !go up until vacancy is found
            do while(cj_attached(endpath)+no_EPS(endpath)>0)
                endpath=endpath+L
                if(endpath>L*H)then
                	write(*,*) "Out of the system, out of luck"
                        stop
                end if
            end do
            if(cj_attached(endpath)+no_EPS(endpath)>0)then
		write(*,*) "endpath wrong"
		stop
	    end if
        !find the distance in xcor and ycor from i to endpath
            signx=1
            signy=1
            xdis=0
            ydis=int(endpath/L)-int(i/L)
	    if(mod(endpath,L)==0 .and. mod(i,L)/=0)then
                ydis=ydis-1
	    else if(mod(endpath,L)/=0 .and. mod(i,L)==0)then
		ydis=ydis+1
	    end if
            if(mod(endpath,L)/=0 .and. mod(i,L)/=0)then
                xdis=mod(endpath,L)-mod(i,L) !this can be negative, or 0
            else if(mod(endpath,L)==0 .and. mod(i,L)/=0)then
                xdis=L-mod(i,L)
            else if(mod(endpath,L)/=0 .and.mod(i,L)==0)then
                xdis=mod(endpath,L)-L
            end if

            if(xdis<0)then
                xdis=abs(xdis)
                signx=-1
            end if

            if(ydis<0)then
                ydis=abs(ydis)
                signy=-1
            end if
            shovnow=0 !first we shove the excreted EPS particle into the site in the path
            do while(xdis>0 .or. ydis>0)
                movex=0
                if(rand()>0.5 .or. ydis==0)then
                    movex=1
                    if(xdis==0)then
                        movex=0
                    end if
                end if
                if(movex==1)then !we shove 1 step in the x direction
                   i=i+signx
                   if(i>L*H)then
                   	write(*,*) "got out of borders x"
                        stop

                   end if
           !see what is already there
                    if(cj_attached(i)+no_EPS(i)==0)then !If we stumble upon a vacancy on the way, just fill the vacancy and end shoving
                        xdis=1
                        ydis=0
                        endpath=i
                    else if(no_EPS(i)>0)then
                        shovnext=0
                        no_EPS(i)=no_EPS(i)-1
                    else if(cj_attached(i)>0)then
                        shovnext=1
                        cj_attached(i)=cj_attached(i)-1
                    end if
!IF there is nothing there, just finish the path early. We have encountered another vacancy on the path, no need to continue shoving.
                    if(shovnow==1)then
                        cj_attached(i)=cj_attached(i)+1
                    else if(shovnow==0)then
                        no_EPS(i)=no_EPS(i)+1
                    end if
                    shovnow=shovnext
                    xdis=xdis-1
                else if(movex==0)then !we shove 1 step in the y direction
                    i=i+signy*L
                   if(i>L*H)then
                   	write(*,*) "got out of borders y"
                        stop

                   end if
                   if(cj_attached(i)+no_EPS(i)==0)then
                        xdis=0
                        ydis=1
                        endpath=i
			
                    else if(no_EPS(i)>0)then
                        shovnext=0
                        no_EPS(i)=no_EPS(i)-1
                    else if(cj_attached(i)>0)then
                        shovnext=1
                        cj_attached(i)=cj_attached(i)-1
                    end if

                    if(shovnow==1)then
                        cj_attached(i)=cj_attached(i)+1
                    else if(shovnow==0)then
                        no_EPS(i)=no_EPS(i)+1
                    end if
                    if(cj_attached(i)+no_EPS(i)>1)then
                        write(*,*) "Unlucky. Tight squeeze", "attached",cj_attached(i),"EPS", no_EPS(i), i
                        stop
                    end if
                    shovnow=shovnext
                    ydis=ydis-1
                end if
            end do
        end if
	rls_invaded(rls,i)=1 !we save the invaded location
	if(Ncj_att/=sum(cj_attached))then
        	write(*,*) "boohoo", Ncj_att, sum(cj_attached), xmin, xmax, i, endpath
                stop
        end if


    end subroutine EPS_generation

!---------------------------OTHER MINOR ROUTINES--------------------------------------------------------
    subroutine Qnext_choosing_ECM_uptake_YES
       gcj=0
        EPS=0
        dbs=0
        dbl=0
        do i=1,L*H
            if(cj_attached(i)>0)then
                if(concentration(i)>0)then
                    if(i<L+1)then
                    	gcj(i)=cj_attached(i)*gmax*(concentration(i)/(concentration(i)+Kc))!Monod equation for growth on the surface
			if(gcj(i)<0)then
				gcj(i)=0
			end if
                    else
                    	gcj(i)=cj_attached(i)*gmax*(concentration(i)/(concentration(i)+Kc)) !growth above the surface
                    end if
                    tot_bonus=0
                    k=1
                    do while (k<5)
                    	if(finalnnodes(i,k)>0)then
                        	if(no_EPS(finalnnodes(i,k))>=(1/NEupt))then
                                        tot_bonus=tot_bonus+bonus
                                        no_EPS(finalnnodes(i,k))=no_EPS(finalnnodes(i,k))-(1/NEupt)
                                end if
                        end if
                        k=k+1
                    end do
                    gcj(i)=gcj(i)+tot_bonus
                    EPS(i)=cj_attached(i)*EPS_max*(1-(concentration(i)/(concentration(i)+Kc)))
                end if
                if(uptk_c(i)<umin_c)then
                    dbs(i)=dbs_cj*cj_attached(i)
                end if
                if((concentration_o(i)>0 .or. concentration(i)>0) .and. cj_attached(i)>0)then
                    oxy=concentration_o(i)
                    conc=concentration(i)
                    dbl(i)=(Al*concentration(i)+Betty*concentration_o(i))*cj_attached(i)

                end if
            end if
        end do
        Q=0
        Qnext=0
        Q1=sum(gcj) !C.jejuni growth
        Q2=sum(gcj)+sum(EPS) !EPS generation
        Q3=sum(gcj)+sum(EPS)+sum(dbs) !Death by starvation
        Q4=sum(gcj)+sum(EPS)+sum(dbs)+sum(dbl) !Death by autolysis
        Q=Q4 !Total rate of events
        

        Qnext=Q*rand()!establishing the next event taking place

    end subroutine Qnext_choosing_ECM_uptake_YES
    subroutine Qnext_choosing
       gcj=0
        EPS=0
        dbs=0
        dbl=0
        do i=1,L*H
            if(cj_attached(i)>0)then
                if(concentration(i)>0)then
                    if(i<L+1)then
                    	gcj(i)=cj_attached(i)*gmax*(concentration(i)/(concentration(i)+Kc))!Monod equation for growth on the surface
						if(gcj(i)<0)then
							gcj(i)=0
						end if
                    else
                    	gcj(i)=cj_attached(i)*gmax*(concentration(i)/(concentration(i)+Kc)) !growth above the surface
                    end if
                    EPS(i)=cj_attached(i)*EPS_max*(1-(concentration(i)/(concentration(i)+Kc)))
                end if
                if(uptk_c(i)<umin_c)then
                    dbs(i)=dbs_cj*cj_attached(i)
                end if
                if((concentration_o(i)>0 .or. concentration(i)>0) .and. cj_attached(i)>0)then
                    oxy=concentration_o(i)
                    conc=concentration(i)
                    dbl(i)=(Al*concentration(i)+Betty*concentration_o(i))*cj_attached(i)

                end if
            end if
        end do
        Q=0
        Qnext=0
        Q1=sum(gcj) !C.jejuni growth
        Q2=sum(gcj)+sum(EPS) !EPS generation
        Q3=sum(gcj)+sum(EPS)+sum(dbs) !Death by starvation
        Q4=sum(gcj)+sum(EPS)+sum(dbs)+sum(dbl) !Death by autolysis
        Q=Q4 !Total rate of events
        

        Qnext=Q*rand()!establishing the next event taking place

    end subroutine Qnext_choosing

    subroutine realtime_increment
        if(Q==0)then
            dt=0
            write(*,*) "nothing will happen ever again", sum(cj_attached)
            stop
        else
            if (t_random==1) then
            	dt=log(1-rand())*(-1)*(1/Q)
            else
				dt=(1/Q)
            end if
        end if
        if(tstep==1)then
            t(tstep)=dt !time increment
        else
            t(tstep)=t(tstep-1)+dt
        end if
        waiting_t(rls,tstep)=dt

    end subroutine realtime_increment

    subroutine death_by_starvation
        i=0
        eventcount=Q2
        do while(eventcount<Qnext .and. i<L*H )
            i=i+1
            eventcount=eventcount+dbs(i) ! dbs_event(i)=cj_attached(i)*dbs_local
        end do
        rls_ev_loc(rls,tstep)=i
        Ncj=Ncj-1
        Ncj_att=Ncj_att-1
        cj_attached(i)=cj_attached(i)-1
        cj_deactivated(i)=cj_deactivated(i)+1


      


    end subroutine death_by_starvation


    subroutine autolysis
        i=0
        eventcount=Q3
        do while(eventcount<Qnext .and. i<L*H )
            i=i+1
            eventcount=eventcount+dbl(i) 
        end do
        rls_ev_loc(rls,tstep)=i
        Ncj=Ncj-1
        Ncj_att=Ncj_att-1
        cj_attached(i)=cj_attached(i)-1
        no_EPS(i)=no_EPS(i)+1

    end subroutine autolysis

    subroutine finding_summit
        if(Ncj_att>0)then
            summit=L*H !starting from the top of the system, scanning downwards for the biofilm peak
            do while(cj_attached(summit)+no_EPS(summit)==0)
                summit=summit-1
            end do
            if(mod(summit,L)==0)then
                summit=int(summit/L)*L+(L)
            else
                summit=int(summit/L)*L+(2*L) !row above the summit
            end if
            if(summit>=(H-1)*L)then
                write(*,*) "boundary reached"
                tstep=maxstep
                write(*,*) "reason 4", summit
            end if
        else
            summit=0
        end if
    end subroutine finding_summit

    subroutine SE_Calc
    no_pores=0   
    do i=1,summit
        if(cj_attached(i)+no_EPS(i)==0)then ! we will only check the empty sites under the summit to see if they are pores
        	next=1
                last=1
                queue=0 !4*L*H array
                checked=0 !L*H array
                queue(1)=i
                top=0
                pore=0
                do while(top==0 .and. pore==0) !this will stop either when we either reach the top or proclaim a site a pore
                	m=1
                        !start of one expanding step
                        do while(next<last+1 .and. top==0) !this will end as long as there are new elements in queue forming, or when an element we are checking is at the summit. Stopping
                           !start of checking one element in the queue
                           	k=1
				do while (k<5 .and. top==0)
                                        if (finalnnodes(k,queue(next))>summit)then
						top=1 !reached the summit
                            		else if(finalnnodes(k,queue(next))>0 .and. finalnnodes(k,queue(next))<summit+1 )then
                                		if(cj_attached(finalnnodes(k,queue(next)))+no_EPS(finalnnodes(k,queue(next)))==0)then
                                                	if(checked(finalnnodes(k,queue(next)))==0)then
                                    				queue(last+m)=finalnnodes(k,queue(next)) !adding a new empty node to check at the end of the queue
                                    				if(queue(last+m)>summit)then
                                        				top=1 !reached the summit, i.e. not a pore
                                   				end if
                                    				m=m+1
							end if
                                		end if
                            		end if
                                        k=k+1
                           	end do
                           !end of checking one element in the queue
                           	checked(queue(next))=1 !making sure we always check the outwards neighbours rather than looking back inwards, to make sure the area of inspection increases in size
                            	next=next+1 !going to check the next element in the queue, until we reach the last element during one step
                        end do
                        !end of one expanding step
                        next=last+1 !at the new step, will go forward in the queue
                        last=next+m-1 ! new last element of the list
                        if(queue(next)==0 .and. top==0)then !if the path is blocked on all sides and there are no new elements in the queue to check, site i is a pore
                           pore=1
                           no_pores=no_pores+1
                        end if
                    end do

                end if
      	end do

    end subroutine SE_Calc

     subroutine nutrient_distribution
        if(sum(cj_attached)==0)then
            concentration=c !in the simplest case, when there is no biofilm, the nutrient concentration will simply stay constant
            concentration_o=o
        else !if there is biofilm in the system, the concentration field must be calculated
!Calculating biofilm density to calculate the mean diffusion
        b_area=0
        biofilm_density=0
        do i=1,summit
        	if(cj_attached(i)+no_EPS(i)>0)then
                        b_area=b_area+1 !this summation does not count pores
			biofilm_density=biofilm_density+cj_attached(i)+no_EPS(i) !making the sum (will divide by b_area after the loop)
		end if
	end do
        !Adding the pores to biofilm area
        b_area=b_area+no_pores
        biofilm_density=biofilm_density/b_area
!Diffusion calculation         
        Diff(tstep)=diff_noeps*(1-biofilm_density)+diff_eps*biofilm_density    
	Alpha=0
	Alpha_o=0
        do i=1,summit
        	if(cj_attached(i)>0)then
                    Alpha(i)=(uptk_c(i)/Diff(tstep))*cj_attached(i)*(dx)**2
                end if
                if(cj_attached(i)>0)then
                   Alpha_o(i)=(uptk_o(i)/Diff(tstep))*cj_attached(i)*(dx)**2
                end if

        end do


            !------Computing the inverse of M-------------------
        !  concentration_matrix=Mat !First assume that the whole field is occupied by biofilm (no boundary layer)
        ! do i=summit+1,L*H  !Now, change the concentration matrix, taking the boundary layer into account
            !    Alpha(i)=c
             !   Alpha_o(i)=o
              !  concentration_matrix(i,:)=0 ! 0 the whole boundary layer row
               ! concentration_matrix(i,i)=1 ! 1 ensures constant concentration in the boundary layer
            !end do
           ! Matrix=concentration_matrix
            !M_inv=0
           ! do d=1,L*H
            !    M_inv(d,d)=1 !Starting with the unity matrix, the inverse matrix will develop
           ! end do
            !do d=1,L*H
             !   dentry=Matrix(d,d) ! dentry is never 0 :) that's the nature of this Matrix
              !  Matrix(d,:)=Matrix(d,:)/dentry  !dividing a row by the diagonal entry
               ! M_inv(d,:)=M_inv(d,:)/dentry !dividing the inverse matrix by the diagonal entry of M
              !  do r=1,L*H
               !     if(r/=d)then
              !         rentry=Matrix(r,d)
               !         if(rentry/=0)then
                !            Matrix(r,:)=Matrix(r,:)+(-1)*rentry*Matrix(d,:) !making the significant row entry 0
                 !           M_inv(r,:)=M_inv(r,:)+(-1)*rentry*M_inv(d,:) !updating the inverse matix
                  !      end if
                   ! end if
                !end do
            !end do
            !-------SIGNIFICANT MATRIX PART INVERSE CALCULATION
	    
            concentration_matrix=0
	    do i=1,summit
                concentration_matrix(i,:)=Mat(i,:)
            end do
            Matrix=concentration_matrix
            M_inv=0
            do d=1,summit
                M_inv(d,d)=1 !Starting with the unity matrix, the inverse matrix will develop
            end do
            do d=1,summit
                dentry=Matrix(d,d) ! dentry is never 0 :) that's the nature of this Matrix
                Matrix(d,:)=Matrix(d,:)/dentry  !dividing a row by the diagonal entry
                M_inv(d,:)=M_inv(d,:)/dentry !dividing the inverse matrix by the diagonal entry of M
                do r=1,summit
                    if(r/=d)then
                        rentry=Matrix(r,d)
                        if(rentry/=0)then
                            Matrix(r,:)=Matrix(r,:)+(-1)*rentry*Matrix(d,:) !making the significant row entry 0
                            M_inv(r,:)=M_inv(r,:)+(-1)*rentry*M_inv(d,:) !updating the inverse matix
                        end if
                    end if
                end do
            end do


!-------End of finding inverse

!Computing the concentration field

            concentration=c ! Concentration field c before the calculations
            concentration_o=o
            do i=1,summit
                concentration(i)=(sum(M_inv(i,:)*Alpha(:)))+c
                if(concentration(i)<0)then
                    concentration(i)=0 ! in case conc is negative
                end if
                concentration_o(i)=(sum(M_inv(i,:)*Alpha_o(:)))+o
                if(concentration_o(i)<0)then
                    concentration_o(i)=0 ! in case conc is negative
                end if
		!introducing scaling up of chemical concentrations due to biofilm porosity
		conc=concentration(i)
		oxy=concentration_o(i)
                uptk_c(i)=UmaxC*concentration(i)/(concentration(i)+Kc)
                uptk_o(i)=UmaxO*concentration_o(i)/(concentration_o(i)+Ko)
            end do
            
        end if


!-------Checking that the calculated inverse is in fact an inverse
        !    multi=0
         !   do row=1,L*H
          !      do col=1,L*H
           !         multi(row,col)=sum(concentration_matrix(row,:)*M_inv(:,col))
            !    end do
          !  end do
           ! do row=1,L*H
            !    if(multi(row,row)>1.01 .or. multi(row,row)<0.99)then
             !      write(*,*) "inverse error 1", multi(row,row), tstep, concentration_matrix(row,row), M_inv(row,row)
              !     stop
              !  end if
               ! do col=1,L*H
                !    if(row/=col)then
                 !       if(multi(row,col)>0.01 .or. multi(row,col)<-0.01)then
                  !          write(*,*) "inverse error 2", multi(row,col)
                   !         stop
              !          end if
               !     end if
               ! end do
           ! end do

!First, found the sites at which the concentration is constant. These are the sites above the biofilm summit (1 site length separation distance

    end subroutine nutrient_distribution

    subroutine Normal_Deviate(gd1,gd2) !subroutine generating random numbers from the Normal distribution
    	real :: zr1,zr2,xr,yr,rho_gauss,gd1,gd2,fac
  

        2101 continue
  	zr1=rand()
  	zr2=rand()
  	xr= 2.0d0*zr1 -1.0d0
  	yr= 2.0d0*zr2 -1.0d0
 	rho_gauss=xr*xr+yr*yr
  	if((rho_gauss.gt.1.0d0).or.(rho_gauss.eq.0.0d0))goto 2101
  
  	fac=sqrt((-2.d0*log(rho_gauss))/rho_gauss)
  
 	gd1=xr*fac
 	gd2=yr*fac
  
 	return
 
    end subroutine Normal_Deviate


end program main
