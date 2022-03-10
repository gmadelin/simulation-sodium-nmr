% =================================================================================================
% function    : prepare_delay_2001
% -------------------------------------------------------------------------------------------------
% purpose     : get date and time for dataname
% input       : -
% output      : date_time_num = numerical date_time string 
% comment     : -  
% reference   : - 
% -------------------------------------------------------------------------------------------------
% date-author : 2015/12 - guillaume.madelin@nyumc.org
%               2018/05 - guillaume.madelin@nyumc.org 
% =================================================================================================
function [date_time_num] = prepare_get_date_time_0201()

    % ---- get date and time
    p.date.date     = datestr(clock);
    p.date.date_vec = datevec(p.date.date);
    
    % ---- numbers -> string
    p.date.year     = num2str(p.date.date_vec(1)-2000);
    p.date.month    = num2str(p.date.date_vec(2));
    p.date.day      = num2str(p.date.date_vec(3));
    p.date.hour     = num2str(p.date.date_vec(4));
    p.date.minute   = num2str(p.date.date_vec(5));
    p.date.second   = num2str(p.date.date_vec(6));
    
    % ---- add '0' in name for numbers < 10 (for uniform formatting)
    if p.date.date_vec(1)<10,  p.date.year   = ['0' p.date.year  ]; end
    if p.date.date_vec(2)<10,  p.date.month  = ['0' p.date.month ]; end
    if p.date.date_vec(3)<10,  p.date.day    = ['0' p.date.day   ]; end
    if p.date.date_vec(4)<10,  p.date.hour   = ['0' p.date.hour  ]; end
    if p.date.date_vec(5)<10,  p.date.minute = ['0' p.date.minute]; end
    if p.date.date_vec(6)<10,  p.date.second = ['0' p.date.second]; end
    
    % ---- final string = YearMonthDay_HourMinuteSecond format
    p.date.date_num = [p.date.year p.date.month p.date.day '_' p.date.hour p.date.minute p.date.second];
    date_time_num = p.date.date_num;
    
end
% =================================================================================================


