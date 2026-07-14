function mesh = mark_elements(adapt, mesh)
% function mark_elements(adapt)
% marks elements for refinement according to 
% the marking strategy stated in 'adapt'
% Possible strategies are
% GR: global (uniform) refinement,  
% MS: maximum strategy,  
% Doerfler: Doerfler 'bulk' strategy

mesh.mark = zeros(mesh.n_elem, 1);

switch adapt.strategy
case 'GR'
  mesh.mark = adapt.n_refine * ones(size(mesh.estimator));
case 'MS'
  mesh.mark(find(mesh.estimator > adapt.MS_gamma*mesh.max_est)) = adapt.n_refine;
case 'Doerfler'
  est_sum2_marked = 0;
  threshold = adapt.Doerfler_theta^2 * mesh.est_sum2;
  gamma = 1;
  while (est_sum2_marked < threshold)
    gamma = gamma - adapt.Doerfler_nu;
    f = find(mesh.estimator > gamma * mesh.max_est);
    mesh.mark(f) = adapt.n_refine;
    est_sum2_marked = sum((mesh.estimator(f)).^2);
  end
end
