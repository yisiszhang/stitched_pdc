function regionTable = ibl_load_region_table(cacheRoot)
%IBL_LOAD_REGION_TABLE Load cached Alyx brain-region metadata from ONE cache.

restDir = fullfile(cacheRoot, '.rest');
assert(exist(restDir, 'dir') == 7, 'Could not find .rest cache under %s', cacheRoot);

files = dir(fullfile(restDir, '*'));
files = files(~[files.isdir]);
[~, ord] = sort([files.bytes], 'ascend');
files = files(ord);

target = '';
for k = 1:numel(files)
    pathk = fullfile(files(k).folder, files(k).name);
    if files(k).bytes > 5e6
        continue;
    end
    txt = fileread(pathk);
    if contains(txt, '"acronym": "void"') && contains(txt, '"parent"') ...
            && contains(txt, '"related_descriptions"') && contains(txt, '"count":')
        target = pathk;
        break;
    end
end

assert(~isempty(target), 'Could not locate cached brain-region metadata in %s', restDir);

blob = jsondecode(fileread(target));
results = blob{1}.results;
n = numel(results);

ids = zeros(n,1);
acronyms = strings(n,1);
parents = nan(n,1);
names = strings(n,1);
for k = 1:n
    ids(k) = results(k).id;
    acronyms(k) = string(results(k).acronym);
    if isempty(results(k).parent)
        parents(k) = NaN;
    else
        parents(k) = results(k).parent;
    end
    names(k) = string(results(k).name);
end

regionTable = table(ids, acronyms, parents, names, ...
    'VariableNames', {'id', 'acronym', 'parent', 'name'});
end
