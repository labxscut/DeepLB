% reads_prob_file has one header line
% Each line is a read. All columns are delimited by TAB. There is one header line.
% Column 1:  marker_index, 1-based index. Those bins that are not markers, are indexed as 0.
% Column 2+: each column is a probability value of a specific tissue. For example, given 16 tissues, there will be 16 columns.
function [theta,q] = em(reads_prob_file, num_tissues, max_iter)
fid = fopen(reads_prob_file);
fmt = ['%d' repmat(' %f',1,num_tissues)];
C = textscan(fid,fmt,'HeaderLines',1);
fclose(fid);

num_reads = length(C{1,1});
p = [];
for i=1:num_tissues,
    p = [p, C{1,i+1}];
end

% initialize q and theta
q = zeros(num_reads, num_tissues);
theta = ones(1,num_tissues) / num_tissues;

for iter=1:max_iter,
% E-step, update q
    q = repmat(theta,num_reads,1).*p;
    q = norm_row(q);
    % M-step, update theta
    theta = sum(q,1);
    theta = norm_row(theta);
    fprintf('theta (round %d): ', iter);
    fprintf('%f  ', theta);
    fprintf('\n');
end

end

function Y = norm_row(X)
row_sums = sum(X,2);
Y = X./repmat(row_sums,1,size(X,2));
end

function Y = norm_col(X)
col_sums = sum(X,1);
Y = X./repmat(col_sums,size(X,1),1);
end

