! cut events requiring the sac files (one-day-long segments)
program main
use sacio
implicit none
type (sac_head) :: sachead1,sachead2,sacheado
integer nn,nstmax,neqmax
parameter (nn=4000000,nstmax=1000,neqmax=5000)
integer out_dd,ieq
integer i,j,nsta,error,year,n1,ic
integer year_b,day_b,year_e,day_e,dhour
integer nlen,nzhour,nzmin,nzsec,nzmsec,nerr,nzero
integer is,ih,iy,id,dayb,daye,jday,icom
integer dsec,month,day,hour,minn,sec,msec,ntemp
integer monthday(12)/31,28,31,30,31,30,31,31,30,31,30,31/
character(8)sta(nstmax)
character(2)net(nstmax)
character(3)com(3),co
character(10)nd,year_day
character(180)command,sac_out
character(80)para,list,eq_list
character(80)sac,dir,sac1,sac2,sac_tmp
character(80)dir_sac,dir_out
real sig(nn),sigall(nn),dt,beg,stla,stlo,sigo(nn)
real evla,evlo,evdp,mag,deltat
logical ext,ext1,ext2

if(iargc().ne.1)then
   write(*,*)'Usage: cut_events param.dat'
   write(*,'("param.dat:")')
   write(*,'("station list")')
   write(*,'("earthquake list")')
   write(*,'("icom component dsec number of components;component;length of output file in sec")')
   write(*,'("dir_of_SAC")')
   write(*,'("dir_of_output")')
   call exit(-1)
endif
i=1
call getarg(1,para)
open(9,file=para)
read(9,'(a80)')list
read(9,'(a80)')eq_list
read(9,*)icom,co,dsec
read(9,'(a80)')dir_sac
read(9,'(a80)')dir_out
!read(9,*)seed_type
write(command,'("mkdir -p",1x,1a,1x,"2>/dev/null")')trim(dir_out)
call system(command)
close(9)
if(icom.eq.1)then
   com(1)=co
else
   com(1)=trim(co)//'Z'
   com(2)=trim(co)//'N'
   com(3)=trim(co)//'E'
endif
open(10,file=list)
do i=1,nstmax
   read(10,*,end=12,err=12)net(i),sta(i)
enddo
12 close(10)
nsta=i-1
open(12,file=eq_list)
do ieq=1,neqmax
   read(12,*,end=14,err=14)year,month,day,hour,minn,sec,msec,evla,evlo,evdp,mag
   write(*,'(7i5.4,1x,4f9.4)')year,month,day,hour,minn,sec,msec,evlo,evla,evdp,mag
   monthday(2)=28
   if(mod(year,4).eq.0.and.mod(year,100).ne.0.or.mod(year,400).eq.0)monthday(2)=29 !leap year
   jday=day
   if(month.gt.1)then
      do i=1,month-1
         jday=jday+monthday(i)
      enddo 
   endif 
   do is=1,nsta
      do ic=1,icom
         write(sac1,'(1a,"/",i4.4,"_",i3.3,"/",i4.4,"_",i3.3,"_00_",1a,"_",1a,"_",1a,".SAC")')&
         trim(dir_sac),year,jday,year,jday,net(is),trim(sta(is)),com(ic)
         write(sac1,'(1a,"/",i4.4,"_",i3.3,"/",i4.4,"_",i3.3,"_00_00_00_",1a,"_",1a,"_",1a,".SAC")')&
         trim(dir_sac),year,jday,year,jday,net(is),trim(sta(is)),com(ic)
         write(sac2,'(1a,"/",i4.4,"_",i3.3,"/",i4.4,"_",i3.3,"_00_",1a,"_",1a,"_",1a,".SAC")')&
         trim(dir_sac),year,jday+1,year,jday+1,net(is),trim(sta(is)),com(ic)
         write(sac_out,'(1a,"/",i4.4,"_",i2.2,"_",i2.2,"_",i2.2,"_",i2.2,"_",i2.2,"_"&
         ,i3.3,"_",1a,"_",1a,".SAC")')trim(dir_out),year,month,day,hour,minn,sec,msec,&
         trim(sta(is)),com(ic)
         write(*,*)trim(sac1),' ',trim(sac2)
         inquire(file=trim(sac1),exist=ext1)
         inquire(file=trim(sac_out),exist=ext2)
         if (.not.ext1.or.ext2)cycle
         write(*,*)'read sac file ',trim(sac1)
         write(*,*)'write to file ',trim(sac_out)
         call read_sachead(sac1,sachead1,nerr)
         if(nerr.eq.-1)cycle
         call read_sac(sac1,sigall,sachead1,nerr)
         out_dd=int(dsec/sachead1%delta)+1
         if(nerr.eq.-1)cycle
         ! deltat=mod(msec,int(dt*1000))/1000.0
         n1=int((hour*60.0*60.0+minn*60.0+sec*1.0+msec*1.0/1000.0)/sachead1%delta)+1
         write(*,*)'n1=',n1
         sig=0
         do i=1,out_dd
            sig(1:out_dd)=sigall(n1:n1+out_dd)
         enddo
         if(n1+out_dd.gt.sachead1%npts)then
            ntemp=n1+out_dd-sachead1%npts
            !call readsac(trim(sac2),sigall,nlen,dt,beg,stla,stlo,nerr) 
            call read_sachead(sac2,sachead2,nerr)
            call read_sac(sac2,sigall,sachead2,nerr)
            if(nerr.ne.-1)sig(sachead1%npts-n1:out_dd)=sigall(1:ntemp)
         endif
         !call taper(sig,out_dd)
         sacheado=sachead1
         sacheado%npts=out_dd
         sacheado%nzyear=year
         sacheado%nzjday=jday
         sacheado%nzmin=minn
         sacheado%nzsec=sec
         sacheado%nzmsec=msec
         sacheado%o=0
         sacheado%mag=mag
         sacheado%evdp=evdp
         sacheado%evlo=evlo
         sacheado%evla=evla
         call write_sac(sac_out,sig,sacheado,nerr)
        enddo ! end loop over component
    enddo  !loop over stations
enddo     ! loop over events
14 close(12)
end program
