function vars = shared_variables(datapath)
    vars.nTrials = 120;
    vars.H = .15;
    vars.replace = true;
    vars.sigma_c = 20;%ob
    vars.sigma_x = 10;%cp
    vars.sigma_y = 33;%28.9516;
    vars.nDegrees = 360;
    vars.path = datapath;
    vars.generate_data = false;
    vars.average = false;
    %"/Users/LabManager/Dropbox (Brown)/arousalLearningPerception/vwm_task/dataForModel/2001_allBlockData.mat";
    save("shared_variables.mat");
end 
