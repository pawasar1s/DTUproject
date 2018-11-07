function mpc = toggle_OID(mpc, on_off)
%TOGGLE_OID Enable, disable or check status of fixed reserve requirements.
% MPC = TOGGLE_OID(MPC, on)
% MPC = TOGGLE_OID(MPC, off)
% T_F = TOGGLE_OID(MPC, status)
if strcmp(upper(on_off), 'ON')
    % <code to check for required OID fields in mpc>
    % add callback functions
    mpc = add_userfcn(mpc, 'ext2int', @userfcn_OID_ext2int);
    mpc = add_userfcn(mpc, 'formulation', @userfcn_OID_formulation);
    mpc = add_userfcn(mpc, 'int2ext', @userfcn_OID_int2ext);
    mpc = add_userfcn(mpc, 'printpf', @userfcn_OID_printpf);
    mpc = add_userfcn(mpc, 'savecase', @userfcn_OID_savecase);
    mpc.userfcn.status.dcline = 1;
elseif strcmp(upper(on_off), 'OFF')
    mpc = remove_userfcn(mpc, 'savecase', @userfcn_OID_savecase);
    mpc = remove_userfcn(mpc, 'printpf', @userfcn_OID_printpf);
    mpc = remove_userfcn(mpc, 'int2ext', @userfcn_OID_int2ext);
    mpc = remove_userfcn(mpc, 'formulation', @userfcn_OID_formulation);
    mpc = remove_userfcn(mpc, 'ext2int', @userfcn_OID_ext2int);
    mpc.userfcn.status.dcline = 0;
elseif strcmp(upper(on_off), 'STATUS')
    if isfield(mpc, 'userfcn') && isfield(mpc.userfcn, 'status') && ...
            isfield(mpc.userfcn.status, 'dcline')
        mpc = mpc.userfcn.status.dcline;
    else
        mpc = 0;
    end
else
    error('toggle_dcline: 2nd argument must be on, off or status');
end