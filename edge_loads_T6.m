function elod = edge_loads_T6(coord, B, eps1)
    L1 = [1;4;1]/6; stride = 2;

    top = find(abs(coord(:,2)-B) < eps1);
    top = sortrows([top, coord(top,1)], 2); top = top(:,1);
    nel = numel(top); ept = [top(1:stride:nel-1); top(nel)];
    lx  = coord(ept,1); dl = lx(2:end) - lx(1:end-1);
    re  = zeros(nel,1);
    for k = 1:numel(dl)
        loc = stride*(k-1) + (1:stride+1);
        re(loc) = re(loc) + dl(k)*L1;
    end
    elod_top = [top, re];

    bot = find(abs(coord(:,2)+B) < eps1);
    bot = sortrows([bot, coord(bot,1)], 2); bot = bot(:,1);
    nel = numel(bot); ept = [bot(1:stride:nel-1); bot(nel)];
    lx  = coord(ept,1); dl = lx(2:end) - lx(1:end-1);
    re  = zeros(nel,1);
    for k = 1:numel(dl)
        loc = stride*(k-1) + (1:stride+1);
        re(loc) = re(loc) + dl(k)*L1;
    end
    elod = [[elod_top(:,1); bot], [elod_top(:,2); -re]];
end
