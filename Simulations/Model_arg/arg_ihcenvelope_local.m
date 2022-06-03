function definput=arg_ihcenvelope_local(definput)
%ARG_IHCENVELOPE
%   #License: GPL
%   #Author: Peter Soendergaard (2011): Initial version
%   #Author: Alejandro Osses (2020): Extensions
%   #Author: Piotr Majdak (2021): Adapted to AMT 1.0
 
  definput.flags.ihc = {'ihc','no_ihc'}; 
  definput.flags.ihc_type={'ihc_undefined','ihc_bernstein1999','ihc_breebaart2001','ihc_dau1996','hilbert', ...
                    'ihc_lindemann1986','ihc_meddis1990','ihc_king2019','ihc_relanoiborra2019','ihc_hwr'};

  definput.keyvals.ihc_minlvl=[];
  definput.keyvals.ihc_filter_order=1;
  definput.keyvals.ihc_scal_constant=[];
  
  definput.groups.ihc_breebaart2001={'ihc_filter_order',5};
  definput.groups.ihc_dau1996={'ihc_filter_order',1};

