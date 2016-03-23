#if defined(_SWE_SIMD)

#include "Compilation_control.f90"

MODULE SWE_SIMD

	use SFC_edge_traversal

	PUBLIC SWE_SIMD_geometry

	TYPE t_SWE_SIMD_geometry
	
		INTEGER :: d ! number of divisions on each edge of the triangular domain (input)
		INTEGER :: n ! number of triangular cells in the patch --> d^2 + 3*d (INCLUDES GHOST CELLS!!)
		! first d*d cells are normal cells. Edges from (d*d+1) to (d*d+3*d) are ghost cells, which are numbered 
		! in counterclockwise order: first d ghost cells belong to the left leg, d+1 to 2*d belong to the hypotenuse,
		! and 2*d+1 to 3d belong to the right leg.
		
		! arrays describing edges. 
		! This means that triangles with ids edges_a[i] and edges_b[i] are neighbors, and edge number i is between them.
		INTEGER, DIMENSION(_SWE_SIMD_NUM_EDGES) :: edges_a, edges_b
		! Describes edge orientation: 1 = parallel to left leg, 2 = to hypotenuse, 3 = to right leg
		INTEGER, DIMENSION(_SWE_SIMD_NUM_EDGES)	:: edges_orientation
		
		! coordinates of cells vertices, used for producing visual output.
		REAL (kind = GRID_SR), DIMENSION(2,3,(_SWE_SIMD_ORDER*_SWE_SIMD_ORDER)) :: coords
		
		! one vector of 2x2 transformation matrices for each i_plotter_type, to avoid repeated computation of this
		REAL (kind = GRID_SR), DIMENSION(_SWE_SIMD_NUM_EDGES_ALIGNMENT,2,2,-8:8) :: transform_matrices

		! represents relationships between cells within a coarse patch and its two refined children
		REAL (kind = GRID_SR), DIMENSION(_SWE_SIMD_ORDER_SQUARE,2) :: first_child, second_child

		contains
		
		procedure, pass :: init => SWE_SIMD_geometry_init
		procedure, pass :: compute_horizontal_internal_edges => SWE_SIMD_compute_horizontal_internal_edges
		procedure, pass :: compute_diagonal_internal_edges => SWE_SIMD_compute_diagonal_internal_edges
		procedure, pass :: compute_ghost_edges => SWE_SIMD_compute_ghost_edges
		procedure, pass :: compute_coords => SWE_SIMD_compute_coords
		procedure, pass :: compute_transform => SWE_SIMD_compute_transform
		procedure, pass :: compute_adaptivity => SWE_SIMD_compute_adaptivity
		procedure, pass	:: cell_in_which_the_point_lies

	END TYPE

	type(t_SWE_SIMD_geometry) SWE_SIMD_geometry

	contains 
	
	!*****************************************************************
	! Procedures/functions for t_SWE_SIMD_geometry
	!*****************************************************************
	
	! procedure SWE_SIMD_geometry%init(d)
	SUBROUTINE SWE_SIMD_geometry_init(geom, d)
		CLASS(t_SWE_SIMD_geometry), INTENT(INOUT) :: geom
		INTEGER, INTENT(IN) :: d
		
		INTEGER :: edges_computed
		INTEGER :: i
		
		! compute number of nodes/triangles (including ghost cells), internal, boundary (ghost) and total edges
		geom%d = d
		geom%n = d*d + 3*d
		
		! compute edges --> they are divided into horizontal, diagonal and ghost edges, but all go into the same arrays
		edges_computed = 0

		call geom%compute_horizontal_internal_edges(edges_computed)
		call geom%compute_diagonal_internal_edges(edges_computed)
		call geom%compute_ghost_edges(edges_computed)
		
		! TODO: (or TOTEST): maybe it would be better to sort the list of edges, might lead to a better memory access pattern
		
		
		! compute cells' vertices' coordinates 
		call geom%compute_coords()
		
		! compute transformation matrices for the solver
		call geom%compute_transform()

		! compute relationships between parent cells and children cells
		call geom%compute_adaptivity()
	
	END SUBROUTINE	

	SUBROUTINE SWE_SIMD_compute_horizontal_internal_edges(geom,edges_computed)
		IMPLICIT NONE
		CLASS(t_SWE_SIMD_geometry), INTENT(INOUT) :: geom
		INTEGER, INTENT(INOUT) :: edges_computed
		INTEGER :: i, last_of_row, step, d
		
		d = geom%d
		last_of_row = 1
		step = 2 !size of step between two neighbor triangles (in different rows)
		i=1
		DO WHILE (i <= (d-1)*(d-1)) 
			edges_computed = edges_computed + 1
			!print *, "edge: ", i, " - ", i+step
			geom%edges_a(edges_computed) = i
			geom%edges_b(edges_computed) = i+step
			geom%edges_orientation(edges_computed) = 2 ! these edges are always parallel to the hypotenuse
			IF (i == last_of_row) THEN
				last_of_row = last_of_row + step + 1
				step = step + 2
				i = i + 1
			ELSE
				i = i + 2
			END IF
		END DO		
	END SUBROUTINE
	
	SUBROUTINE SWE_SIMD_compute_diagonal_internal_edges(geom,edges_computed)
		IMPLICIT NONE
		CLASS(t_SWE_SIMD_geometry), INTENT(INOUT) :: geom
		INTEGER, INTENT (INOUT) :: edges_computed
		INTEGER :: i, edge_type
		
		edge_type = 3 ! these edges' orientation will be 1 or 3, alternatedly
 		DO i=1,geom%d * geom%d
 			IF (isPerfectSquare(i) == 0) THEN
				edges_computed = edges_computed + 1
				geom%edges_a(edges_computed) = i
				geom%edges_b(edges_computed) = i+1
				geom%edges_orientation(edges_computed) = edge_type
				edge_type = mod(edge_type+2, 4) ! 1 becomes 3, 3 becomes 1
 			END IF
 		END DO
		
		
	END SUBROUTINE

	!returns true if n is a perfect square, false otherwise
	integer function isPerfectSquare(n)
		IMPLICIT NONE
		integer , intent (IN) :: n
		real :: root

		!if the square root has no decimal portion
		!test this with 'if sqrt == truncated decimal sqrt'
		root = sqrt(real(n))
		if (root == aint(root)) then
			isPerfectSquare = 1
		else
			isPerfectSquare = 0
		end if
	end function
	
	SUBROUTINE SWE_SIMD_compute_ghost_edges(geom,edges_computed)
		IMPLICIT NONE
		CLASS(t_SWE_SIMD_geometry), INTENT(INOUT) :: geom
		INTEGER, INTENT (INOUT) :: edges_computed
		INTEGER :: i, j, cell, d
		
		d = geom%d
		
		! left boundary --> 1st of i_th row is neighbor of cell d*d + i
		DO i=1,d ! i means row here
			cell = (i-1)*(i-1)+1
			edges_computed = edges_computed + 1
			geom%edges_a(edges_computed) = d*d+i
			geom%edges_b(edges_computed) = cell
			geom%edges_orientation(edges_computed) = 1
 		END DO
 		
 		! hypotenuse (bottom) boundary
 		i = (d-1)*(d-1)+1 ! 1st cell of last row
 		j = d*d+d+1 ! 1st cell of hypotenuse ghost cells
 		DO WHILE (i<=d*d)
			edges_computed = edges_computed + 1
			geom%edges_a(edges_computed) = i
			geom%edges_b(edges_computed) = j
			geom%edges_orientation(edges_computed) = 2
 		
			i = i + 2
			j = j + 1
 		END DO
 		
		! right boundary --> last of i_th row is neighbor of cell d*d + 3*d + 1 - i
		DO i=1,d ! i means row here
			cell = i*i
			edges_computed = edges_computed + 1
			geom%edges_a(edges_computed) = cell
			geom%edges_b(edges_computed) = d*d + 3*d + 1 - i
			geom%edges_orientation(edges_computed) = 3
 		END DO
		
		
	END SUBROUTINE
	
	! computes the coordinates of the vertices in a local [0,1]x[0,1] representation. 
	! The (patch) triangle legs lie on the x and y axes. 
	! The output traversals typically apply transformations on these coordinates.
	SUBROUTINE SWE_SIMD_compute_coords(geom)
		IMPLICIT NONE
		CLASS(t_SWE_SIMD_geometry), INTENT(INOUT) :: geom

		integer											:: i, row, col
		real (kind = GRID_SR), dimension(2,3)			:: r_points
		real (kind = GRID_SR), dimension(2)				:: tmp
		real (kind = GRID_SR)							:: d
		
		! length of each leg edge
		d = 1.0_GRID_SR/_SWE_SIMD_ORDER
		row=1
		col=1
		
		! first cell
		r_points(:,1) = [0.0_GRID_SR, 	0.0_GRID_SR	]
		r_points(:,2) = [0.0_GRID_SR, 	d			]
		r_points(:,3) = [d, 			0.0_GRID_SR	]
		
		do i=1,_SWE_SIMD_ORDER * _SWE_SIMD_ORDER
			
			! copy coordinates to array
			geom%coords(:,:,i) = r_points
			
			! make sure that all triangles are described in counterclockwise order
			if (mod(col,2) == 1) then
				tmp = geom%coords(:,1,i)
				geom%coords(:,1,i) = geom%coords(:,2,i)
				geom%coords(:,2,i) = tmp
			end if

			!compute next points
			col = col + 1
			if (col == 2*row) then
				col = 1
				row = row + 1
				r_points(:,1) = [d*(row-1), 	0.0_GRID_SR	]
				r_points(:,2) = [d*(row-1), 	d	 		]
				r_points(:,3) = [ d*row, 		0.0_GRID_SR ]
			else
				r_points(:,3) = r_points(:,1)
				r_points(:,1) = r_points(:,2)
				if (mod(col,2) == 0) then
					r_points(1,2) = r_points(1,2) - d
				else
					r_points(2,2) = r_points(2,2) + d
				end if
			end if
		end do
	
	END SUBROUTINE
	
	SUBROUTINE SWE_SIMD_compute_transform(geom)
		IMPLICIT NONE
		CLASS(t_SWE_SIMD_geometry), INTENT(INOUT) :: geom

		integer											:: i, j
		real(kind = GRID_SR)							:: normals(2,3)
		
		do i = -8, 8 ! for each i_plotter_type
            if (i .ne. 0) then
				! copy/compute normal vectors
				! normal for type 2 edges is equal to the 2nd edge's normal
				normals(:,2) = ref_plotter_data(i)%edges(2)%normal 
				! normal for types 1 and 3 depend on cell orientation.
				! notice that normal for type 1 points inwards. That's why it is multiplied by -1.
				if (i < 0) then ! orientation = backward
					normals(:,1) = - ref_plotter_data(i)%edges(1)%normal 
					normals(:,3) = ref_plotter_data(i)%edges(3)%normal 
				else ! orientation = forward, so reverse edge order
					normals(:,1) = - ref_plotter_data(i)%edges(3)%normal 
					normals(:,3) = ref_plotter_data(i)%edges(1)%normal 				
				end if
				
				! compute transformations matrices
				do j=1,_SWE_SIMD_NUM_EDGES 
					geom%transform_matrices(j,1,:,i) = normals(:,geom%edges_orientation(j))
					geom%transform_matrices(j,2,:,i) = [ - normals(2,geom%edges_orientation(j)), normals(1,geom%edges_orientation(j)) ]
				end do
            end if
        end do
	END SUBROUTINE
	
    ! finds out relationships between cells on refinement and coarsening operations
	SUBROUTINE SWE_SIMD_compute_adaptivity(geom)
		IMPLICIT NONE
		CLASS(t_SWE_SIMD_geometry), INTENT(INOUT) 		:: geom
		integer											:: i, j
		real (kind = GRID_SR)							:: d, square_of_sqrt2_div_2
		real (kind = GRID_SR), DIMENSION(2,2)			:: transform_1, transform_2
		real (kind = GRID_SR), DIMENSION(2)				:: center_point, test1, test2

		! length of each leg edge
		d = 1.0_GRID_SR/_SWE_SIMD_ORDER

		! used for computing transformations
		square_of_sqrt2_div_2 = sqrt(2.0_GRID_SR)*0.5_GRID_SR
		square_of_sqrt2_div_2 = square_of_sqrt2_div_2 * square_of_sqrt2_div_2

		! considering that the original patch has vertices (0,0), (1,0), (0,1),
		! the first child will be (0,0), (1,0), (sqrt(2),sqrt(2))
		! and the second child will be (0,0), (sqrt(2),sqrt(2)), (0,1)

		! === compute for first child ===
		! transformation = rotate by (-135) degrees and scale by sqrt(2)/2
		! sin(-135) = -sqrt(2)/2, cos(-135) = -sqrt(2)/2
		transform_1(1,:) = [-square_of_sqrt2_div_2, +square_of_sqrt2_div_2]
		transform_1(2,:) = [-square_of_sqrt2_div_2, -square_of_sqrt2_div_2]

		! === compute for second child ===
		! transformation = rotate by (+135) degrees and scale by sqrt(2)/2
		! sin(+135) = sqrt(2)/2, cos(+135) = -sqrt(2)/2
		transform_2(1,:) = [-square_of_sqrt2_div_2, -square_of_sqrt2_div_2]
		transform_2(2,:) = [+square_of_sqrt2_div_2, -square_of_sqrt2_div_2]

		do i=1, _SWE_SIMD_ORDER_SQUARE ! i = cell id in child
			! compute baricenter of cell
			center_point = [0.0_GRID_SR, 0.0_GRID_SR]
			do j=1,3
				center_point = center_point + geom%coords(:,j,i)
			end do
			center_point = center_point/3.0_GRID_SR

			! a child cell will tipically lie with one or two cells belonging to the parent
			! --> test a point slightly to the left and other slightly to the right and store in which
			! 	  parent cell they lie. The result can be the same for both, but it doesn't matter.
			test1 = center_point + [0.1_GRID_SR*d, 0.0_GRID_SR]
			test2 = center_point + [-0.1_GRID_SR*d, 0.0_GRID_SR]
			! apply transformation --> a translation and then the computed transform matrix
			test1 = test1 - [1.0_GRID_SR, 0.0_GRID_SR]
			test2 = test2 - [1.0_GRID_SR, 0.0_GRID_SR]
			test1 = matmul(transform_1,test1)
			test2 = matmul(transform_1,test2)

			! save the relationships -> parent cells in first_child(i,:) intersect first child cell #i and so on
			geom%first_child(i,1) = geom%cell_in_which_the_point_lies(test1)
			geom%first_child(i,2) = geom%cell_in_which_the_point_lies(test2)

			! do the same for the second child --> note that the translation is also different
			test1 = center_point + [0.1_GRID_SR*d, 0.0_GRID_SR]
			test2 = center_point + [-0.1_GRID_SR*d, 0.0_GRID_SR]
			test1 = test1 - [0.0_GRID_SR, 1.0_GRID_SR]
			test2 = test2 - [0.0_GRID_SR, 1.0_GRID_SR]
			test1 = matmul(transform_2,test1)
			test2 = matmul(transform_2,test2)

			geom%second_child(i,1) = geom%cell_in_which_the_point_lies(test1)
			geom%second_child(i,2) = geom%cell_in_which_the_point_lies(test2)
		end do
	END SUBROUTINE

	integer FUNCTION cell_in_which_the_point_lies(geom,point)
		CLASS(t_SWE_SIMD_geometry), INTENT(INOUT) 				:: geom
		real (kind=GRID_SR), dimension(2), intent(in)			:: point
		integer													:: i

		do i=1, _SWE_SIMD_ORDER_SQUARE
			if (is_point_inside_triangle(point, geom%coords(:,1,i), geom%coords(:,2,i), geom%coords(:,3,i)) > 0) then
				cell_in_which_the_point_lies = i
				exit
			end if
		end do
	end FUNCTION

	! returns whether a->b->c is a left turn
	integer FUNCTION left_turn(a, b, c)
		real (kind=GRID_SR), dimension(2), intent(in)			:: a, b, c
		
		if ( (b(1) - a(1)) * (c(2) - b(2)) - (b(2) - a(2)) * (c(1) - b(1)) > 0 ) then
			left_turn = 1
		else
			left_turn = 0
		end if
		
	end FUNCTION
	
	! returns whether the given point is inside the COUNTER-CLOCKWISE triangle defined by the given vertices a, b, c
	integer FUNCTION is_point_inside_triangle(point, a, b, c)
		IMPLICIT NONE
		real (kind=GRID_SR), dimension(2), intent(in)			:: point, a, b, c		
		
		if ( left_turn(a,b,point)>0 .and. left_turn(b,c,point)>0 .and. left_turn(c,a,point)>0 ) then
			is_point_inside_triangle = 1
		else
			is_point_inside_triangle = 0
		end if

	END FUNCTION
	
END MODULE
#endif
