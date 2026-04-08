addpath(genpath(pwd));

theta1_list_deg = [-30, -40, -50];

% Dugdale fit constants chosen for each initial angle
CD_list = [17.2, 1.3*17.2, 2.0*17.2];

AllResults = cell(numel(theta1_list_deg),1);

for i = 1:numel(theta1_list_deg)
    C0 = cfg_nonperforated();

    C0.study.eta_list = 80:5:120;
    C0.study.save_results = false;
    C0.run.verbose = true;

    % initial crack angle
    C0.geometry.theta1 = deg2rad(theta1_list_deg(i));

    % Dugdale-based delta sizing
    C0.study.use_dugdale_delta = true;
    C0.study.dugdale_C = CD_list(i);

    Results = run_eta_kinking_study(C0);

    AllResults{i} = Results;
end