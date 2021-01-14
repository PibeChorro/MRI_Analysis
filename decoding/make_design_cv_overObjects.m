

function design = make_design_cv_overObjects(cfg)

%% generate design matrix (CV)
design.function.name = mfilename;
design.function.ver = 'v2020_new';


design.set = [1 1 1];

tricks                  = ceil(cfg.files.chunk/2);
train_cards_sticks      = tricks~=1;
test_balls              = tricks==1;

train_balls_sticks      = tricks~=2;
test_cards              = tricks==2;

train_balls_cards       = tricks~=3;
test_sticks             = tricks==3;

design.train        = double(horzcat(train_cards_sticks,   train_balls_sticks, train_balls_cards));
design.test         = double(horzcat(test_balls,           test_cards,         test_sticks));
design.label        = double(repmat (cfg.files.label,1,3));

end