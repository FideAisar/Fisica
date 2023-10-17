
c      PARAMETER(mxnphs=20)  !maximum number of available phase spaces.
      PARAMETER(mxnphs=76)  !maximum number of available phase spaces.
                         ! Those actually implemented are enumerated 
                         ! at the beginning of proc.f 

*common for single process
* igroup=group number of the selected process according to the table 
*        in amp.f
* n_kphs=total number of phase spaces active in the selected process 
* iphs_ind = phase space index which identifies a selected  channel 
*             in the process at hand 
* ialfa(i)=0/1 off/on phase space i for the selected process 
* idp(j)= particle identity according to the amplitude order 
*         for the selected process
* idp_inout(j)=-1/+1  incoming/outgoing particles according 
*              to the amplitude order for the selected process
* iorder(j,i) gives the momenta order as passed to the amplitude
*             from the phase space i of the selected process
*             The order of the particles is different for the 
*             different phase spaces.
*             iorder(j,i)=k means that the k-th particle in the
*             order of the amplitude is the j-th in the order
*             of the i-th phase space               
* identicalinitial=1/2  number of identical particle in the initial 
*                       state



      integer iphs_ind,ialfa(mxnphs),
     &        idp(8),idp_inout(8),iorder(8,mxnphs),nidenticalinitial,
     &        igo(8,6),igr(6),ngluon,ngin

      COMMON/fxninput_proc/n_kphs,iphs_ind,ialfa,
     &      idp,idp_inout,iorder,nidenticalinitial,igo,igr,ngluon,ngin

*sandro 25/07
* common for the process as read from readinput
      COMMON/process/iproc(8)
*sandro 25/07 end

* common for the integration:        
*  rnormal(i) normalization value for a given channel
*  alfa(i)  value of the alfas of a multichannel integration for 
*  a given process 
* init and it are flags of vegas which are used in fxn.f integ.f
*   and oneshot.f

      COMMON/subinteg/rnormal(mxnphs),alfa(mxnphs)
      COMMON/flat_integ/init,it
 
