function [date_itrf_changes] = readITRFChanges(path)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

raw_text = fileread(path);
text_cell = strsplit(raw_text, '\n');
date_itrf_changes = datetime(text_cell, 'InputFormat', 'yyyy-MM-dd', 'TimeZone', 'UTC');

end

