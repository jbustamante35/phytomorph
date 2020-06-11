#!/bin/bash
##
echo 'Adding Collaborators'
collaborator='tina.miller'
echo $collaborator
./addCollaborator "$collaborator"
./listCollaborator
./addCollaborator 'george miller'
./listCollaborator
##
echo 'Adding Services'
service1='hello.world'
./addService $service1
./listService
./addService 'second.service'
./listService
##
echo 'Adding XchangSite'
xchangeName='point1'
./addXchangesite $xchangeName $collaborator
./listXchange
##
echo 'Linking site to service'
./addServiceToPoint $xchangeName $service1
