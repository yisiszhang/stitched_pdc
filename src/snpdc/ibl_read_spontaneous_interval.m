function interval = ibl_read_spontaneous_interval(csvFile)
%IBL_READ_SPONTANEOUS_INTERVAL Read spontaneous interval from IBL passive CSV.

raw = readcell(csvFile, 'Delimiter', ',');
headers = string(raw(1, :));
col = find(headers == "spontaneousActivity", 1);
if isempty(col)
    interval = [];
    return;
end

rowNames = string(raw(2:end, 1));
startRow = find(rowNames == "start", 1) + 1;
stopRow = find(rowNames == "stop", 1) + 1;

if isempty(startRow) || isempty(stopRow)
    interval = [];
    return;
end

interval = [raw{startRow, col}, raw{stopRow, col}];
interval = double(interval);
end
