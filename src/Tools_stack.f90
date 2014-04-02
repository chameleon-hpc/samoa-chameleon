! Sam(oa)² - SFCs and Adaptive Meshes for Oceanic And Other Applications
! Copyright (C) 2010 Oliver Meister, Kaveh Rahnema
! This program is licensed under the GPL, for details see the file LICENSE


!> Generic container template
!> Warning: this a template module body and requires preprocessor commands before inclusion.
!> Usage: Inside your module definition, insert
!>
!> #define _CNT_DATA_TYPE			<value>
!> #define _CNT_TYPE_NAME			<value>
!>
!> #include <this_file>
!>
!> where
!>
!> @item _CNT_DATA_TYPE				Any derived type, base data type for the container
!> @item _CNT_CHUNK_SIZE			Any positive integer, sets the capacity of the vector chunks for vector-type stacks
!> @item _CNT_TYPE_NAME				Container type name, will be used prefix for all operations
!>
!> The resulting container is defined as <container_name>, the chunk as <container_name>_chunk, methods as <container_name>_<method>
!> @author Oliver Meister

!multiple levels of indirection are necessary to properly resolve the names
#define _CONC2(X, Y)			X ## _ ## Y
#define _PREFIX(P, X)			_CONC2(P, X)
#define _CNT_(X)				_PREFIX(_CNT_TYPE_NAME, X)

#define _CNT					_CNT_TYPE_NAME
#define _T						_CNT_DATA_TYPE

PRIVATE create, destroy, is_empty, push_element, push_pointer, pop_element, pop_pointer, current_pointer, to_string
!PRIVATE push_elements, pop_elements

!define a common interface for all derived types

type _CNT
	_T, dimension(:), pointer	:: elements => null()		!< element array
	integer (kind = GRID_SI)	:: i_current_element = 0	!< number of elements stored on stack
    logical                     :: exists = .true.

	contains

	procedure, pass :: create
	procedure, pass :: destroy
	procedure, pass :: is_empty
	procedure, pass :: reset

	procedure, pass :: push_data => push_element
	procedure, pass :: pop_data => pop_element
	!procedure, pass :: push_elements
	!procedure, pass :: pop_elements

	procedure, pass :: push => push_pointer
	procedure, pass :: current => current_pointer
	procedure, pass :: pop => pop_pointer

	procedure, pass :: to_string
end type _CNT

contains

!> Creates a new stack
subroutine create(stack, i_capacity)
	class(_CNT), intent(inout)					:: stack					!< stack object
	integer (kind = GRID_SI), intent(in)		:: i_capacity				!< stack capacity

	integer (kind = GRID_SI)					:: i_error

    assert(.not. associated(stack%elements))
	allocate(stack%elements(i_capacity), stat = i_error); assert_eq(i_error, 0)

	stack%i_current_element = 0
end subroutine create

!> Destroys the stack
subroutine destroy(stack)
	class(_CNT), intent(inout)					:: stack					!< stack object

	integer (kind = GRID_SI)					:: i_error

	assert(associated(stack%elements))
	deallocate(stack%elements, stat = i_error); assert_eq(i_error, 0)
end subroutine

!array access

subroutine push_element(stack, data)
	class(_CNT), intent(inout)					:: stack					!< stack object
	_T, intent(in)								:: data						!< data

	stack%i_current_element = stack%i_current_element + 1

	!check for stack overflow
	assert_ge(stack%i_current_element, lbound(stack%elements, 1))
	assert_le(stack%i_current_element, ubound(stack%elements, 1))

	stack%elements(stack%i_current_element) = data
end subroutine

!subroutine push_elements(stack, data)
!	class(_CNT), intent(inout)					:: stack					!< stack object
!	_T, intent(in)								:: data(:)					!< data
!
!	stack%i_current_element = stack%i_current_element + size(data)
!
!	!check for stack overflow
!	assert_ge(stack%i_current_element, lbound(stack%elements, 1))
!	assert_le(stack%i_current_element, ubound(stack%elements, 1))
!
!	stack%elements(stack%i_current_element - size(data) + 1 : stack%i_current_element) = data
!end subroutine

subroutine pop_element(stack, data)
	class(_CNT), intent(inout)					:: stack					!< stack object
	_T, intent(out)								:: data						!< data

	!check for stack overflow
	assert_ge(stack%i_current_element, lbound(stack%elements, 1))
	assert_le(stack%i_current_element, ubound(stack%elements, 1))

	data = stack%elements(stack%i_current_element)

	stack%i_current_element = stack%i_current_element - 1
end subroutine

!subroutine pop_elements(stack, data)
!	class(_CNT), intent(inout)					:: stack					!< stack object
!	_T, intent(out)								:: data(:)					!< data
!
!	!check for stack overflow
!	assert_ge(stack%i_current_element, lbound(stack%elements, 1))
!	assert_le(stack%i_current_element, ubound(stack%elements, 1))
!
!	data = stack%elements(stack%i_current_element : stack%i_current_element - size(data) + 1 : -1)
!
!	stack%i_current_element = stack%i_current_element - size(data)
!end subroutine

function push_pointer(stack) result(p_data)
	class(_CNT), intent(inout)					:: stack					!< stack object
	_T, pointer									:: p_data					!< data pointer

	stack%i_current_element = stack%i_current_element + 1

	!check for stack overflow
	assert_ge(stack%i_current_element, lbound(stack%elements, 1))
	assert_le(stack%i_current_element, ubound(stack%elements, 1))

	p_data => stack%elements(stack%i_current_element)
end function

function current_pointer(stack) result(p_data)
	class(_CNT), intent(inout)					:: stack					!< stack object
	_T, pointer									:: p_data					!< data pointer

	!check for stack overflow
	assert_ge(stack%i_current_element, lbound(stack%elements, 1))
	assert_le(stack%i_current_element, ubound(stack%elements, 1))

	p_data => stack%elements(stack%i_current_element)
end function

function pop_pointer(stack) result(p_data)
	class(_CNT), intent(inout)					:: stack					!< stack object
	_T, pointer									:: p_data					!< data pointer

	!check for stack overflow
	assert_ge(stack%i_current_element, lbound(stack%elements, 1))
	assert_le(stack%i_current_element, ubound(stack%elements, 1))

	p_data => stack%elements(stack%i_current_element)

	stack%i_current_element = stack%i_current_element - 1
end function

!> checks if the stack is empty
function is_empty(stack)
	class(_CNT), intent(in)					    :: stack					!< stack object
	logical										:: is_empty

	assert_ge(stack%i_current_element, 0)

	is_empty = (stack%i_current_element == 0)
end function

elemental function to_string(stack) result(str)
	class(_CNT), intent(in)						:: stack					!< stack object
	character (len = 32)						:: str

	if (size(stack%elements) > 0) then
		write(str, '(A, I0, A, I0)') "capacity: ", size(stack%elements), " current: ", stack%i_current_element
	else
		write(str, '(A, I0)') "capacity: ", size(stack%elements)
	endif
end function

elemental subroutine reset(stack)
	class(_CNT), intent(inout)					    :: stack					!< stack object

	stack%i_current_element = 0
end subroutine


#undef _CNT_DATA_TYPE
#undef _CNT_TYPE_NAME
