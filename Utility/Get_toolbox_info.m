function info_toolbox = Get_toolbox_info(script_str)
% function info_toolbox = Get_toolbox_info(script_str)

[ver,ver_str] = fastACI_version;
info_toolbox.version = ver;
info_toolbox.version_str = ver_str;
info_toolbox.date = Get_date_and_time_str;
info_toolbox.mfilename = script_str;
