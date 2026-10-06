import ViaB

/-! Print transitive axiomatic dependencies of every public theorem.
Expected: only Lean/mathlib's standard propext, Classical.choice, Quot.sound;
no sorryAx and no custom geometric axioms. Explicit hypotheses remain visible
in the declarations and are separately discussed in README.md. -/
#print axioms ViaB.curvaturePolynomial_pos
#print axioms ViaB.screenGap_pos
#print axioms ViaB.directionB_at_horizon
#print axioms ViaB.directionB_at_infinity_negative
#print axioms ViaB.wronskian_hasDerivAt
#print axioms ViaB.ratio_hasDerivAt
#print axioms ViaB.wronskian_strictMonoOn
#print axioms ViaB.ratio_deriv_pos
#print axioms ViaB.derivative_nonpos_at_first_zero
#print axioms ViaB.no_first_zero
#print axioms ViaB.positive_right_of_positive_derivative
#print axioms ViaB.exists_first_zero_of_nonpos
#print axioms ViaB.sturm_positive
#print axioms ViaB.sturm_positive_at_reference_zero
#print axioms ViaB.trace_tangent_nondegenerate
#print axioms ViaB.covering_fiberCount_locallyConstant
#print axioms ViaB.covering_injective_of_one_fiber
#print axioms ViaB.curvaturePolynomial_negative_below_threshold
#print axioms ViaB.screenGap_negative_below_threshold
#print axioms ViaB.screenGap_radial
#print axioms ViaB.directionB_continuousAt
#print axioms ViaB.directionB_le_at_one
#print axioms ViaB.directionB_uniform_far
#print axioms ViaB.directionB_large_radius
#print axioms ViaB.radialSpeed_at_infinity_pos
#print axioms ViaB.radialSpeed_le_charge
#print axioms ViaB.radialSpeed_uniform_far
#print axioms ViaB.increment_ge_of_derivative_lower
#print axioms ViaB.increment_le_of_derivative_upper
#print axioms ViaB.positive_right_of_positive_value
#print axioms ViaB.positive_pair_invariant
#print axioms ViaB.excursion_log_bound
#print axioms ViaB.outgoing_cone_margins
#print axioms ViaB.outgoing_cone_invariant
#print axioms ViaB.outgoing_cone_linear_growth
#print axioms ViaB.radial_growth_tendsto_atTop
#print axioms ViaB.outgoing_cone_escape
#print axioms ViaB.no_zero_limit_of_positive_nondecreasing
#print axioms ViaB.no_radial_horizon_limit
#print axioms ViaB.positive_future_incompatible_with_negative_derivative
#print axioms ViaB.no_incoming_direction_limit
#print axioms ViaB.outward_direction_lt_pi_div_two
#print axioms ViaB.large_excursion_enters_angular_cone
#print axioms ViaB.exists_last_exit
#print axioms ViaB.derivative_nonneg_at_right_exit
#print axioms ViaB.strict_cone_entry
#print axioms ViaB.cone_escape_from_nonstrict_radius
#print axioms ViaB.exists_outgoing_cone_angle
#print axioms ViaB.unbounded_radius_escapes
#print axioms ViaB.bounded_nondecreasing_primitive_converges
#print axioms ViaB.decay_from_convergent_primitive
#print axioms ViaB.decay_from_bounded_primitive
#print axioms ViaB.directionBound_pos
#print axioms ViaB.directionB_abs_bound
#print axioms ViaB.bounded_exterior_finite_angle_sine_decay
#print axioms ViaB.sine_decay_forces_radial_limit
#print axioms ViaB.bounded_exterior_finite_angle_outgoing_limit
#print axioms ViaB.angular_speed_lower
#print axioms ViaB.radialSpeed_lower_near_outgoing
#print axioms ViaB.bounded_radius_converges_from_outgoing_limit
#print axioms ViaB.bounded_monotone_future_converges
#print axioms ViaB.derivative_limit_zero_on_bounded_positive_radius
#print axioms ViaB.zero_outgoing_speed_forces_charge
#print axioms ViaB.no_bounded_exterior_finite_angle
#print axioms ViaB.bounded_exterior_flow_angle_diverges

#check ViaB.sturm_positive
#check ViaB.sturm_positive_at_reference_zero
#check ViaB.covering_injective_of_one_fiber
#check ViaB.outgoing_cone_escape
#check ViaB.large_excursion_enters_angular_cone
#check ViaB.no_radial_horizon_limit
#check ViaB.no_incoming_direction_limit
#check ViaB.unbounded_radius_escapes
#check ViaB.bounded_exterior_flow_angle_diverges

#print axioms ViaB.unbounded_radius_strict_cone_entry
#print axioms ViaB.cone_entry_stable_in_family
#print axioms ViaB.escaping_set_isOpen
#print axioms ViaB.directionField_contDiffAt
#print axioms ViaB.directionField_locally_lipschitz
#print axioms ViaB.trajectory_gap_bound
#print axioms ViaB.trajectory_eval_continuousAt
#print axioms ViaB.directionFlow_state_derivative
#print axioms ViaB.directionFlow_eval_continuousAt
#print axioms ViaB.cone_escape_stable_from_initial_data
#print axioms ViaB.finite_angle_exterior_escapes
#print axioms ViaB.exterior_escape_or_angle_diverges
#print axioms ViaB.trajectory_stays_in_tube
#print axioms ViaB.trajectory_eval_continuousAt_of_local_gap
#print axioms ViaB.cone_escape_stable_of_local_field_bound

#check ViaB.cone_escape_stable_of_local_field_bound
#check ViaB.finite_angle_exterior_escapes

#print axioms ViaB.finite_positive_lower_bound
#print axioms ViaB.compact_uniform_local_lipschitz
#print axioms ViaB.finite_reference_uniform_field_bound
#print axioms ViaB.directionFlow_uniform_field_bound
#print axioms ViaB.directionFlow_finite_time_dependence
#print axioms ViaB.directionFlow_eval_continuousAt_from_initial
#print axioms ViaB.cone_escape_stable_from_regular_reference
#print axioms ViaB.escaping_set_isOpen_from_initial

#check ViaB.finite_reference_uniform_field_bound
#check ViaB.cone_escape_stable_from_regular_reference
#check ViaB.escaping_set_isOpen_from_initial

-- MetricHamiltonian
#print axioms ViaB.raise_lower_identity
#print axioms ViaB.lower_raise_identity
#print axioms ViaB.metric_of_raised_momentum
#print axioms ViaB.null_velocity_ellipse
#print axioms ViaB.cometric_complete_square
#print axioms ViaB.nullAmplitude_pos
#print axioms ViaB.nullAmplitude_sq
#print axioms ViaB.future_null_sheet
#print axioms ViaB.future_null_sheet_unique

-- NullHamiltonianDynamics
#print axioms ViaB.hamiltonian_time_momentum_derivative
#print axioms ViaB.hamiltonian_radial_momentum_derivative
#print axioms ViaB.hamiltonian_angular_momentum_derivative
#print axioms ViaB.momentum_direction_circle
#print axioms ViaB.momentum_direction_sin_pos
#print axioms ViaB.null_radial_velocity
#print axioms ViaB.null_angular_velocity
#print axioms ViaB.vaidyaW_pos
#print axioms ViaB.exterior_iff_vaidyaW_lt_charge
#print axioms ViaB.receiver_timelike
#print axioms ViaB.time_covector_timelike
#print axioms ViaB.hamiltonian_radius_derivative
#print axioms ViaB.radial_momentum_future_reduction
#print axioms ViaB.radial_momentum_direction_combination

-- DirectionFromHamiltonian
#print axioms ViaB.directionCos_hasDerivAt
#print axioms ViaB.directionCos_dynamics
#print axioms ViaB.angle_derivative_from_cosine
#print axioms ViaB.momentumAngle_properties
#print axioms ViaB.momentumAngle_dynamics

-- AffineReparametrization
#print axioms ViaB.affine_hamiltonian_to_direction

#check ViaB.affine_hamiltonian_to_direction
#check ViaB.momentumAngle_dynamics

-- LocalTimeInverse
#print axioms ViaB.local_time_inverse_exists
#print axioms ViaB.affineTimeSpeed_continuousAt
#print axioms ViaB.hamiltonian_time_inverse_exists
#print axioms ViaB.hamiltonian_local_direction_exists

-- TimeDependentHamiltonian
#print axioms ViaB.mass_along_affine_derivative
#print axioms ViaB.hamiltonian_mass_derivative
#print axioms ViaB.hamiltonian_time_coordinate_derivative
#print axioms ViaB.hamiltonian_along_path_derivative
#print axioms ViaB.hamiltonian_flow_energy_derivative_zero
#print axioms ViaB.hamiltonian_null_level_propagates
#print axioms ViaB.hamiltonian_angular_coordinate_derivative
#print axioms ViaB.angular_momentum_constant
#print axioms ViaB.hamiltonian_flow_null_propagates

#check ViaB.hamiltonian_local_direction_exists
#check ViaB.hamiltonian_flow_null_propagates

-- LocalDirectionExistence
#print axioms ViaB.liftedDirectionField_contDiffAt
#print axioms ViaB.direction_local_solution_exists
#print axioms ViaB.direction_local_solution_unique
#print axioms ViaB.direction_exterior_local_solution_exists
#print axioms ViaB.direction_local_initial_value_problem

-- FiniteEndpointLimit
#print axioms ViaB.planar_lipschitz_extension
#print axioms ViaB.bounded_speed_endpoint_limit
#print axioms ViaB.compact_direction_speed_bound
#print axioms ViaB.compact_direction_endpoint_restart

-- RegularEndpointContinuation
#print axioms ViaB.direction_endpoint_left_derivative
#print axioms ViaB.direction_endpoint_backward_agreement
#print axioms ViaB.compact_direction_continues

#check ViaB.direction_local_initial_value_problem
#check ViaB.direction_exterior_local_solution_exists
#check ViaB.compact_direction_continues

-- ExteriorFiniteTube
#print axioms ViaB.direction_radial_derivative
#print axioms ViaB.exterior_w_bounds
#print axioms ViaB.exterior_radial_speed_le_charge
#print axioms ViaB.exterior_finite_radius_upper
#print axioms ViaB.exterior_finite_state_bounds
#print axioms ViaB.exterior_finite_regular_compact_tube
#print axioms ViaB.exterior_finite_coordinate_continuation

-- AngularFiniteInvariance
#print axioms ViaB.direction_angular_derivative
#print axioms ViaB.linear_ode_initial_nonzero
#print axioms ViaB.sine_direction_coefficient_bound
#print axioms ViaB.finite_direction_angle_stays_nonradial
#print axioms ViaB.exterior_finite_angle_invariance
#print axioms ViaB.exterior_nonradial_finite_coordinate_continuation
#print axioms ViaB.exterior_nonradial_endpoint_control
#print axioms ViaB.nonradial_capture_gap_derivative
#print axioms ViaB.exterior_finite_endpoint_capture_or_continues

#check ViaB.exterior_nonradial_endpoint_control
#check ViaB.exterior_finite_endpoint_capture_or_continues
#check ViaB.nonradial_capture_gap_derivative

#print axioms ViaB.direction_continuation_preserves_launch
#print axioms ViaB.exterior_nonradial_initial_value_continuation
#check ViaB.exterior_nonradial_initial_value_continuation


-- IntervalUniqueness
#print axioms ViaB.direction_solution_unique_on_open_connected
#print axioms ViaB.direction_solution_unique_on_interval

-- MaximalExteriorConstruction
#print axioms ViaB.exterior_segments_agree
#print axioms ViaB.maximalExteriorDomain_isOpen
#print axioms ViaB.maximalExteriorDomain_isPreconnected
#print axioms ViaB.maximalExteriorCurve_agrees
#print axioms ViaB.maximalExteriorCurve_solves
#print axioms ViaB.maximalExteriorCurve_preserves_launch
#print axioms ViaB.exterior_segment_exists
#print axioms ViaB.maximalExteriorDomain_eq_Ioo
#print axioms ViaB.maximalExteriorDomain_eq_Ioi
#print axioms ViaB.exterior_segment_domain_subset_maximal
#print axioms ViaB.exterior_global_solution_of_unbounded_endpoints

-- MaximalExteriorCapture
#print axioms ViaB.maximal_exterior_finite_endpoint_is_capture

-- MaximalConeEscape
#print axioms ViaB.outgoing_cone_invariant_on_finite_interval
#print axioms ViaB.maximal_exterior_endpoints_unbounded_after_cone_entry
#print axioms ViaB.maximal_directionFlow_of_unbounded_endpoints
#print axioms ViaB.directionFlow_exists_and_escapes_after_cone_entry
#print axioms ViaB.directionFlow_exists_from_strict_cone_launch

-- DirectionMomentumReconstruction
#print axioms ViaB.reconstructed_nullAmplitude
#print axioms ViaB.reconstructed_direction
#print axioms ViaB.reconstructed_future_null
#print axioms ViaB.reconstructed_time_speed_identity
#print axioms ViaB.reconstructed_radial_momentum_hasDerivAt
#print axioms ViaB.reconstructed_radial_momentum_dynamics
#print axioms ViaB.reconstructed_normalized_hamiltonian_identities

-- ReconstructedTimeMomentum
#print axioms ViaB.time_momentum_force_from_null_radial_equations
#print axioms ViaB.reconstructedPv_differentiableAt
#print axioms ViaB.reconstructed_time_momentum_dynamics

-- ConstructedAffineClock
#print axioms ViaB.continuous_interval_primitive_exists
#print axioms ViaB.local_inverse_derivative_on_neighborhood
#print axioms ViaB.positive_clock_with_inverse_exists
#print axioms ViaB.directionAffineClockRate_pos
#print axioms ViaB.directionAffineClockRate_continuousAt
#print axioms ViaB.direction_affine_clock_exists

-- DirectionHamiltonianLift
#print axioms ViaB.direction_angular_primitive_exists
#print axioms ViaB.direction_to_affine_hamiltonian_at
#print axioms ViaB.direction_affine_hamiltonian_germ_exists
#print axioms ViaB.maximal_exterior_affine_hamiltonian_germ_exists
