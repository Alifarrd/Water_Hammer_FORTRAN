program FSI_MOCQH2
implicit none
INTEGER::i,ne,nt,nx,tt,tx,j,ii
REAL::k,vp,R,E,ee,r_f,r_t,ks,cf,ct,g,cfh,cth,a,b,c,a2,b2,c2,l,dt,zt,dxz,zz,pres,v0,bs,as,hres,q0,x1,x2,az,at,ge
!---------------------------------------------------[guess]
INTEGER::ng,ig,jg,jg1,ig3,jg3,kg
REAL::landag,sg
real,allocatable::ag(:,:),xg(:),bg(:)
!---------------------------------------------------[guess]         
real,allocatable::t(:,:),f(:),gx(:),vi(:,:),pi(:,:),ui(:,:),zi(:,:),bc(:,:)
!OPEN(UNIT=10,FILE='input.txt',STATUS='old',ACTION='read')
!OPEN(UNIT=12,FILE='ans.txt',STATUS='replace',ACTION='write')
OPEN(UNIT=20,FILE='ansplot.plt',STATUS='replace',ACTION='write')
!----------
ng=4
ALLOCATE (ag(ng,ng+1),xg(ng),bg(ng+1))
!----------

!=============================================================[Initial Information]
ge=9.806           !gravitational acceleration(m/s2)
k=2.1e9            !Fluid's Bulk modulus,[pa]
vp=0.29 !0.5!0.29          !Possion ratio,
R=0.797/2          !inner Radius of pipe,[m]
E=2e11             !Young modulus of pipe wall material,[pa]
ee=0.008           !pipe wall thickness,[m]
r_f=1000           !fluid mass density,[kg/m3]
r_t=7870.          !Tube mass density,[kg/m3]
pres=1.e5          !reservoir pressure
Hres=pres/(r_f*ge)
v0=1.002           !initial velosity
!=============================================================
ks=1/((1/k)+((1-vp**2)*((2*r)/(E*ee))))                          !effictive bulk modulus,[pa]
cf=(ks/r_f)**0.5                                                 !Classical Fluid wave speed,[m/s]
ct=(e/r_t)**0.5                                                  !Classical Tube wave speed,[m/s]
g=((1+2*vp**2*(r_f/r_t)*(r/e))*cf**2+ct**2)**0.5                  !constant,[m/s]
cfh=0.5*sqrt(2.)*(g**2-(g**4-4*cf**2*ct**2)**0.5)**0.5            !FSI fluid wave speed,[m/s]
cth=0.5*sqrt(2.)*(g**2+(g**4-4*cf**2*ct**2)**0.5)**0.5            !FSI Tube wave speed,[m/s]
!=============================================================
a=1/(r_f*cfh)
b=2*vp*((cfh/ct)**2/1-(cfh/ct)**2)
c=(2*vp/r_t*cfh)*((cfh/ct)**2/1-(cfh/ct)**2)
!=============================================================
a2=1/(r_t*cth)
b2=((-vp*r*r_f)/(e*r_t))*((cf/cth)**2/1-(cf/cth)**2)
c2=((vp*r)/(e*r_t*cth))*((cf/cth)**2/1-(cf/cth)**2)
!=============================================================
as=(3.14/4)*(2*r)**2
at=(3.14/4)*(2*(r-ee))**2
az=as-at
bs=(cfh)/(ge*as)
q0=as*v0
!=============================================================
l=20                                !length
ne=20                               !number of element
tt=5                                !Total time
nx=ne+1                             !number of node
dxz=l/ne                            
dt=dxz/cfh                         
nt=floor(tt/dt)+1                            !number of time steps
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
ALLOCATE (t(4,5),f(4),gx(4),vi(nx,nt),pi(nx,nt),ui(nx,nt),zi(nx,nt),bc(2,3))
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
t(1,1)=bs
t(1,2)=1
t(1,3)=b*(cfh/ge)
t(1,4)=-c*(cfh/ge)
!----------
t(2,1)=-bs
t(2,2)=1
t(2,3)=-b*(cfh/ge)
t(2,4)=-c*(cfh/ge)
!----------
t(3,1)=-b2/as
t(3,2)=-c2*r_f*ge
t(3,3)=1
t(3,4)=-a2
!----------
t(4,1)=b2/as
t(4,2)=c2*(r_f*ge)
t(4,3)=1
t(4,4)=a2
!======================================================================================[initial boundry]
vi(:,:)=0
pi(:,:)=0
pi(1,:)=hres
ui(:,:)=0
zi(:,:)=0

vi(:,1)=q0
pi(:,1)=hres
ui(:,1)=0
zi(:,1)=hres*(as/az)
!======================================================================================[initial boundry]

do tx=1,nt-1
  
!=====================================================================================[boundry Condition]
  
!===================================[reservoir]
x1=(-bs*vi(2,tx))+(hres-pi(2,tx))-(b*(cfh/ge)*ui(2,tx))+(c*(cfh/ge)*zi(2,tx))
x2=((b2/as)*vi(6,tx))-(c2*r_f*ge*(hres-pi(6,tx)))+(ui(6,tx))+(a2*zi(6,tx))

zi(1,tx+1)=(x2+((b2*x1)/(as*bs)))/(((-b2*c)/(as*bs))+a2)

!zi(1,tx+1)=((pres/(r_f*ge))*(as/az))
vi(1,tx+1)=(x1+(c*zi(1,tx+1)))/(+bs)

!===================================[valve]
       ii=nx;

       pi(ii,tx+1)=pi(ii-1,tx)+(bs*vi(ii-1,tx))-(bs*vi(ii,tx+1))
       zi(ii,tx+1)=(r_f*ge)*(hres-(pi(ii,tx+1)))*(as/az)
       !zi(ii,tx+1)=((r_f*ge)*(pi(ii,tx+1)-pres))*(as/az)
       !zi(ii,tx+1)=(pi(ii,tx+1))*(as/az)
!=====================================================================================[boundry Condition] 
do i=2,nx-1
if(i<=5)then             !....................................................................................... FIRST NODEs:
  zt=dt-(((i-1)*dxz)/cth)
  zz=dt-zt
!-----------------------------------------[opproximation of V,P,U,Z]   
  gx(1)=((zz*vi(1,tx))+(zt*vi(1,tx+1)))/dt 
  gx(2)=((zz*pi(1,tx))+(zt*pi(1,tx+1)))/dt
  gx(3)=((zz*ui(1,tx))+(zt*ui(1,tx+1)))/dt
  gx(4)=((zz*zi(1,tx))+(zt*zi(1,tx+1)))/dt
!-----------------------------------------[opproximation of V,P,U,Z]
 
t(3,5)=gx(3)-(a2*gx(4))-(b2*(cfh/9.81)*gx(1))-(c2*r_f*ge*gx(2))
t(4,5)=gx(3)+(a2*gx(4))+(b2*(cfh/9.81)*gx(1))+(c2*r_f*ge*gx(2))

else if(i>nx-5)then      !....................................................................................... END NODES:
  zt=dt-(((nx-(i-1))*dxz)/cth)
  zz=dt-zt
!-----------------------------------------[opproximation of V,P,U,Z]   
  gx(1)=((zz*vi(nx,tx))+(zt*vi(nx,tx+1)))/dt 
  gx(2)=((zz*pi(nx,tx))+(zt*pi(nx,tx+1)))/dt
  gx(3)=((zz*ui(nx,tx))+(zt*ui(nx,tx+1)))/dt
  gx(4)=((zz*zi(nx,tx))+(zt*zi(nx,tx+1)))/dt
!-----------------------------------------[opproximation of V,P,U,Z]
 
t(3,5)=gx(3)-(a2*gx(4))-(b2*(cfh/ge)*gx(1))-(c2*r_f*ge*gx(2))
t(4,5)=gx(3)+(a2*gx(4))+(b2*(cfh/ge)*gx(1))+(c2*r_f*ge*gx(2))

else                     !.......................................................................................MIDDEL NODES
  
t(3,5)=ui(i-5,tx)-(a2*zi(i-5,tx))-(b2*(cfh/ge)*vi(i-5,tx))-(c2*r_f*ge*pi(i-5,tx))
t(4,5)=ui(i+5,tx)+(a2*zi(i+5,tx))+(b2*(cfh/ge)*vi(i+5,tx))+(c2*r_f*ge*pi(i+5,tx))

end if
  
t(1,5)=pi(i-1,tx)+(bs*vi(i-1,tx))+(b*(cfh/ge)*ui(i-1,tx))-(c*(cfh/ge)*zi(i-1,tx))
t(2,5)=pi(i+1,tx)-(bs*vi(i+1,tx))-(b*(cfh/ge)*ui(i+1,tx))-(c*(cfh/ge)*zi(i+1,tx))
!--------------------------------------------------[guess]
do kg=1,ng
ag(kg,:)=t(kg,:)
end do
do ig=1,ng
!if (ag(ig,ig)==0) then
!do jg=ig+1,ng
!if (ag(jg,ig)/=0) then
!bg=ag(jg,:)
!ag(jg,:)=ag(ig,:)
!ag(ig,:)=bg
!end if
!end do
!end if
do jg1=ig+1,ng
landag=-ag(jg1,ig)/ag(ig,ig)
ag(jg1,:)=(landag*ag(ig,:))+ag(jg1,:)
end do
end do
 xg(ng)=ag(ng,ng+1)/ag(ng,ng)
 do ig3=ng-1,1,-1
 sg=0.
 do jg3=ig3+1,ng
 sg=sg+(ag(ig3,jg3)*xg(jg3))
 end do
 xg(ig3)=(ag(ig3,ng+1)-sg)/ag(ig3,ig3)
 end do
  vi(i,tx+1)=xg(1)
  Pi(i,tx+1)=xg(2)
  ui(i,tx+1)=xg(3)
  zi(i,tx+1)=xg(4)
!--------------------------------------------------[guess]
end do

end do
!---------------------------------------------------[print output]
write(20,*)'titele'
write(20,*)'variables=x,y,vi,hi,ui,zi'
write(20,*)'zone'
write(20,*)'i=',nt,'j=',nx
write(20,*)'f=point'
    do i=1,nx
    do j=1,nt
 write(20,*)i,j,vi(i,j),pi(i,j),ui(i,j),zi(i,j)
 end do
 end do
!---------------------------------------------------
  end program