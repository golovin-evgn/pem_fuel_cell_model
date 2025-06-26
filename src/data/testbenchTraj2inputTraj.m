function u_traj_inp = testbenchTraj2inputTraj(u_traj_tb, p)
% testbenchTraj2inputTraj(u_traj_tb, p) transforms the interpolation
% structure from testbench representation to input representation

data_structure = p.testbench2struct(u_traj_tb.data.');

transformed_data_structure = testbench2model(data_structure,p);

u_traj_inp.time = u_traj_tb.time;
u_traj_inp.data = p.inputs2vec(transformed_data_structure).';

end
