
MODULE SWE_SIMD

	use SFC_edge_traversal

	PUBLIC SWE_SIMD_geometry

	TYPE t_SWE_SIMD_geometry
	
		INTEGER :: d ! number of divisions on each edge of the triangular domain (input)
		INTEGER :: n ! number of triangular cells in the patch --> d^2 + 3*d (INCLUDES GHOST CELLS!!)
		! first d*d cells are normal cells. Edges from (d*d+1) to (d*d+3*d) are ghost cells, which are numbered 
		! in counterclockwise order: first d ghost cells belong to the left leg, d+1 to 2*d belong to the hypotenuse,
		! and 2*d+1 to 3d belong to the right leg.
	
		INTEGER :: num_internal_edges ! number of internal edges --> f(d) is given by the recurrence: f(2) = 3, f(d+1) = f(d) + 3*d
		INTEGER :: num_boundary_edges ! number of edges on the boundary --> d*3 (ghost edges)
		INTEGER :: num_edges ! number of edges (sum of both above)
		
		! arrays describing edges. This means that triangles with ids edges_a[i] and edges_b[i] are neighbors, and edge number i is between them.
		INTEGER, DIMENSION(:), POINTER :: edges_a, edges_b
		
		! coordinates of cells vertices, used for producing visual output.
		REAL (kind = GRID_SR), DIMENSION(2,3,_SWE_SIMD_ORDER*_SWE_SIMD_ORDER) :: coords
		

		contains
		
		procedure, pass :: init => SWE_SIMD_geometry_init
		procedure, pass :: compute_horizontal_internal_edges => SWE_SIMD_compute_horizontal_internal_edges
		procedure, pass :: compute_diagonal_internal_edges => SWE_SIMD_compute_diagonal_internal_edges
		procedure, pass :: compute_ghost_edges => SWE_SIMD_compute_ghost_edges
		procedure, pass :: compute_coords => SWE_SIMD_compute_coords

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
		geom%num_boundary_edges = 3*d
		
		geom%num_internal_edges=3
		DO i=3,d
			geom%num_internal_edges = geom%num_internal_edges +  3*(i-1)
		END DO
		
		geom%num_edges = geom%num_internal_edges + geom%num_boundary_edges
		
	
		! allocate arrays for list of edges
		allocate(geom%edges_a(geom%num_edges))
		allocate(geom%edges_b(geom%num_edges))

		
		! compute edges --> they are divided into horizontal, diagonal and ghost edges, but all go into the same arrays
		edges_computed = 0
		call geom%compute_horizontal_internal_edges(edges_computed)
		call geom%compute_diagonal_internal_edges(edges_computed)
		call geom%compute_ghost_edges(edges_computed)
		
		! TODO: (or TOTEST): maybe it would be better to sort the list of edges, might lead to a better memory access pattern
		
		
		! compute cells' vertices' coordinates 
		call geom%compute_coords()
	
	
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
		INTEGER :: i
		
 		DO i=1,geom%d * geom%d
 			IF (isPerfectSquare(i) == 0) THEN
				edges_computed = edges_computed + 1
				geom%edges_a(edges_computed) = i
				geom%edges_b(edges_computed) = i+1
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
			geom%edges_a(edges_computed) = cell
			geom%edges_b(edges_computed) = d*d+i
 		END DO
 		
 		! hypotenuse (bottom) boundary
 		i = (d-1)*(d-1)+1 ! 1st cell of last row
 		j = d*d+d+1 ! 1st cell of hypotenuse ghost cells
 		DO WHILE (i<=d*d)
			edges_computed = edges_computed + 1
			geom%edges_a(edges_computed) = i
			geom%edges_b(edges_computed) = j
 		
			i = i + 2
			j = j + 1
 		END DO
 		
		! right boundary --> last of i_th row is neighbor of cell d*d + 3*d + 1 - i
		DO i=1,d ! i means row here
			cell = i*i
			edges_computed = edges_computed + 1
			geom%edges_a(edges_computed) = cell
			geom%edges_b(edges_computed) = d*d + 3*d + 1 - i
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
		real (kind = GRID_SR)							:: d
		
		! length of each leg edge
		d = 1.0_GRID_SR/_SWE_SIMD_ORDER
		row=1
		col=1
		
		! first cell
		r_points(1,1) = 0.0
		r_points(2,1) = 0.0
		r_points(1,2) = 0.0
		r_points(2,2) = d
		r_points(1,3) = d
		r_points(2,3) = 0.0
		
		do i=1,_SWE_SIMD_ORDER * _SWE_SIMD_ORDER
			
			! copy coordinates to array
			geom%coords(1,1,i) = r_points(1,1)
			geom%coords(2,1,i) = r_points(2,1)
			geom%coords(1,2,i) = r_points(1,2)
			geom%coords(2,2,i) = r_points(2,2)
			geom%coords(1,3,i) = r_points(1,3)
			geom%coords(2,3,i) = r_points(2,3)
			
			print *, ".", i, geom%coords(:,:,i)
			
			!compute next points
			col = col + 1
			if (col == 2*row) then
				col = 1
				row = row + 1
				r_points(1,1) = d*(row-1)
				r_points(2,1) = 0.0
				r_points(1,2) = d*(row-1)
				r_points(2,2) = d
				r_points(1,3) = d*row
				r_points(2,3) = 0.0
			else
				r_points(1,3) = r_points(1,1)
				r_points(2,3) = r_points(2,1)
				r_points(1,1) = r_points(1,2)
				r_points(2,1) = r_points(2,2)
				
				if (mod(col,2) == 0) then
					r_points(1,2) = r_points(1,2) - d
				else
					r_points(2,2) = r_points(2,2) + d
				end if
			end if
		end do
	
	END SUBROUTINE
	
END MODULE
