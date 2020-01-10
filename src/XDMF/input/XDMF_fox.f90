#include "Compilation_control.f90"
#include "XDMF/XDMF_compilation_control.f90"

! This module provides wrapper functions for commonly used FoX library function
! Essentially, it acts as a namespace to prevent name collisions
module XDMF_fox

    use FoX_dom, only: &
        fox_parseFile               => parseFile, &
        t_fox_Node                  => Node, &
        fox_getDocumentElement      => getDocumentElement, &
        fox_getFirstChild           => getFirstChild, &
        fox_getChildNodes           => getChildNodes, &
        t_fox_NodeList              => NodeList, &
        t_fox_NamedNodeMap          => NamedNodeMap, &
        fox_destroy                 => destroy, &
        fox_getLength               => getLength, &
        fox_getAttributes           => getAttributes, &
        fox_item                    => item, &
        fox_getLocalName            => getLocalName, &
        fox_getNamedItem            => getNamedItem, &
        fox_getNodeValue            => getNodeValue, &
        fox_getLastChild            => getLastChild, &
        fox_extractDataContent      => extractDataContent, &
        fox_extractDataAttribute    => extractDataAttribute, &
        fox_getElementsByTagName    => getElementsByTagName

    contains

end module