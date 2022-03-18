!------------------------------------------------------------------------------
      SUBROUTINE NetcdfOutputSolver_Init0(Model,Solver,dt,Transient)
!------------------------------------------------------------------------------
      USE DefUtils
      IMPLICIT NONE
!------------------------------------------------------------------------------
      TYPE(Solver_t) :: Solver
      TYPE(Model_t) :: Model
      REAL(KIND=dp) :: dt
      LOGICAL :: Transient
!------------------------------------------------------------------------------
      TYPE(ValueList_t), POINTER :: SolverParams

        SolverParams => GetSolverParams()

        CALL ListAddNewLogical( SolverParams,'Optimize Bandwidth',.FALSE.)
        ! current limitations; bulk and no halo
        CALL ListAddNewLogical( SolverParams,'Save Bulk Only',.TRUE.)
        CALL ListAddNewLogical( SolverParams,'Skip Halo Elements',.TRUE.)

!------------------------------------------------------------------------------
      END SUBROUTINE NetcdfOutputSolver_Init0
!------------------------------------------------------------------------------
      
      
!------------------------------------------------------------------------------
      SUBROUTINE NetcdfOutputSolver( Model,Solver,dt,TransientSimulation )
!------------------------------------------------------------------------------
      USE DefUtils
      USE SolverUtils
      USE SaveUtils
      USE Netcdf

#ifdef HAVE_PROJ
      USE fortranc
      USE proj
#endif

      IMPLICIT NONE
      TYPE(Solver_t) :: Solver
      TYPE(Model_t) :: Model
      REAL(dp) :: dt
      LOGICAL :: TransientSimulation

      CHARACTER(*), PARAMETER :: Caller = 'NetcdfOutputSolver'

      INTEGER, SAVE :: nTime = 0
      CHARACTER(MAX_NAME_LEN),SAVE :: NcFile
      INTEGER,SAVE :: ElemFirst, ElemLast
      INTEGER, ALLOCATABLE, TARGET,SAVE :: NodePerm(:),InvNodePerm(:), InvDgPerm(:), DgPerm(:)
      LOGICAL, ALLOCATABLE,SAVE :: ActiveElem(:)
      INTEGER,SAVE :: NumberOfGeomNodes,NumberOfDofNodes,NumberOfElements
      LOGICAL,SAVE :: NoPermutation


      TYPE(Mesh_t), POINTER :: Mesh
      TYPE(ValueList_t),POINTER :: Params
      LOGICAL :: GotIt
      LOGICAL :: Parallel
      INTEGER :: GroupId
      INTEGER :: PEs,Part 

      LOGICAL :: SaveLinear=.TRUE. ! save only corners
      INTEGER :: MaxElementNodes  ! // MaxNumNodesPerFace
      INTEGER,PARAMETER :: Connect_Fill=-1

      CHARACTER(MAX_NAME_LEN) :: BaseFile,OutputDirectory

      LOGICAL :: ComputeLonLat=.FALSE.
#ifdef HAVE_PROJ
      TYPE(pjuv_object) :: coordp,coordg
      TYPE(pj_object) :: pj
      CHARACTER(LEN=MAX_NAME_LEN) :: proj_string
#endif

      Params => GetSolverParams()
      Mesh => Model % Mesh
      ! Comment: Can we get the Mas number of corners instead of nodes?
      MaxElementNodes = Model % MaxElementNodes

      IF (Mesh % MeshDim.NE.2) &
        CALL FATAL(Caller,"Mesh dim should be 2 for now...")

      PEs = ParEnv % PEs
      Part = ParEnv % MyPE
      Parallel = (PEs > 1)

      nTime = nTime + 1
      IF ((nTime>1).AND.(Mesh%Changed)) &
        CALL FATAL(Caller,"mesh has changed; not supported")

      !------------------------------------------------------------------------------
      ! Create and initialise file at first visit
      !------------------------------------------------------------------------------
      IF ( nTime == 1 ) THEN

        !------------------------------------------------------------------------------
        ! Initialize File Naming
        !------------------------------------------------------------------------------
        BaseFile = GetString( Params,'Output File Name',GotIt )
        IF ( .NOT.GotIt ) BaseFile = "NetcdfOutput"

        CALL Info(Caller,'Saving results in Netcdf format with prefix: '//TRIM(BaseFile))
        CALL Info(Caller, 'Saving number of partitions: '//TRIM(I2S(PEs)))

        CALL SolverOutputDirectory( Solver, BaseFile, OutputDirectory, UseMeshDir = .TRUE.  )
        BaseFile = TRIM(OutputDirectory)// '/' //TRIM(BaseFile)

        IF (Parallel) THEN
          IF ( PEs < 10) THEN
           WRITE( NcFile,'(A,A,I1.1,A,I1.1)') TRIM((BaseFile)),"_",PEs,"np",Part
          ELSE IF ( PEs < 100) THEN
           WRITE( NcFile,'(A,A,I2.2,A,I2.2)') TRIM((BaseFile)),"_",PEs,"np",Part
          ELSE IF ( PEs < 1000) THEN
           WRITE( NcFile,'(A,A,I3.3,A,I3.3)') TRIM((BaseFile)),"_",PEs,"np",Part
          ELSE
           WRITE( NcFile,'(A,A,I4.4,A,I4.4)') TRIM((BaseFile)),"_",PEs,"np",Part
          END IF
        ELSE
          NcFile = TRIM((BaseFile))
        END IF
        NcFile = TRIM((NcFile)) // ".nc"

        !------------------------------------------------------------------------------
        ! Initialize stuff for masked saving
        !------------------------------------------------------------------------------
        GroupId = 0
        CALL GenerateSaveMask(Mesh,Params,Parallel,GroupId,SaveLinear,&
               NodePerm,ActiveElem,NumberOfGeomNodes,NumberOfElements,&
               ElemFirst,ElemLast)
        !------------------------------------------------------------------------------
        ! If we have a discontinuous mesh then create the permutation vectors to deal
        ! with the discontinuities.
        !------------------------------------------------------------------------------
        CALL GenerateSavePermutation(Mesh,.FALSE.,.FALSE.,SaveLinear,ActiveElem,NumberOfGeomNodes,&
               NoPermutation,NumberOfDofNodes,DgPerm,InvDgPerm,NodePerm,InvNodePerm)

        ! The partition is active for saving if there are any nodes
        ! to write. There can be no elements nor dofs without nodes.
        CALL ParallelActive( NumberOfDofNodes > 0 )

        IF( ElemLast > Mesh % NumberOfBulkElements ) &
          CALL FATAL(Caller,"Saving boundary elements not supported yet....")

#ifdef HAVE_PROJ
        ComputeLonLat=ListGetLogical(Params,"Compute LonLat")
        IF (ComputeLonLat) THEN
          proj_string = ListGetString(Params,'projection',UnFoundFatal=.TRUE.)
          pj = pj_init_plus(TRIM(proj_string)//CHAR(0))
        END IF
#endif

        CALL CreatNetcdfFile(NcFile,NumberOfDofNodes,NumberOfElements)  
        CALL Info(Caller, 'CreatNetcdfFile Done')

        CALL WriteMeshInfo(NcFile)
        CALL Info(Caller, 'WriteMeshInfo Done')

      END IF

      CALL WriteVariables(NcFile,nTime)

      CONTAINS 
        SUBROUTINE CreatNetcdfFile(FName,NNodes,NBulkElements)
        implicit none
        Character(LEN=MAX_NAME_LEN),INTENT(IN) :: FName
        INTEGER, INTENT(IN) :: NNodes,NBulkElements

        TYPE(Variable_t),POINTER :: Solution
        CHARACTER(LEN=1024) :: Txt, ScalarFieldName
        INTEGER :: Vari,VarType
        INTEGER :: ncid
        INTEGER, DIMENSION(4) :: dimid
        INTEGER :: varid

        integer :: status
        integer :: id

        LOGICAL :: ScalarsExist

        ! do not overwrite file... for now
        call check( nf90_create(TRIM(FName),NF90_NOCLOBBER,ncid))

        call check( nf90_def_dim(ncid,'nMesh_node',NNodes,dimid(1)))
        call check( nf90_def_dim(ncid,'nMesh_face',NBulkElements,dimid(2)))
        call check( nf90_def_dim(ncid,'nMaxMesh_face_nodes',MaxElementNodes,dimid(3)))
        call check( nf90_def_dim(ncid,'time',NF90_UNLIMITED,dimid(4)))

        call check( nf90_def_var(ncid,'Mesh',NF90_INT,varid))
          call check( nf90_put_att(ncid,varid,"cf_role", "mesh_topology") )
          call check( nf90_put_att(ncid,varid,"topology_dimension", 2) )
          IF (ComputeLonLat) THEN
            call check( nf90_put_att(ncid,varid,"node_coordinates", "Mesh_node_lon Mesh_node_lat") )
            call check( nf90_put_att(ncid,varid,"face_coordinates", "Mesh_face_lon Mesh_face_lat") )
          ELSE
            call check( nf90_put_att(ncid,varid,"node_coordinates", "Mesh_node_x Mesh_node_y") )
          ENDIF
          call check( nf90_put_att(ncid,varid,"face_node_connectivity", "Mesh_face_nodes"))
          call check( nf90_put_att(ncid,varid,"face_dimension","nMesh_face"))

        call check( nf90_def_var(ncid,'Mesh_face_nodes',NF90_INT,(/dimid(3),dimid(2)/),varid))
          call check( nf90_put_att(ncid,varid,"cf_role", "face_node_connectivity"))
          call check( nf90_def_var_fill(ncid,varid,0,Connect_Fill))
          call check( nf90_put_att(ncid,varid,"start_index",1))

        call check( nf90_def_var(ncid,'Mesh_node_x',NF90_DOUBLE,dimid(1),varid))
        call check( nf90_def_var(ncid,'Mesh_node_y',NF90_DOUBLE,dimid(1),varid))

        IF (ComputeLonLat) THEN
          call check( nf90_def_var(ncid,'Mesh_node_lon',NF90_DOUBLE,dimid(1),varid))
            call check( nf90_put_att(ncid,varid,"standard_name","longitude"))
            call check( nf90_put_att(ncid,varid,"long_name","nodes longitude"))
            call check( nf90_put_att(ncid,varid,"units","degrees_east"))

          call check( nf90_def_var(ncid,'Mesh_node_lat',NF90_DOUBLE,dimid(1),varid))
            call check( nf90_put_att(ncid,varid,"standard_name","latitude"))
            call check( nf90_put_att(ncid,varid,"long_name","nodes latitude"))
            call check( nf90_put_att(ncid,varid,"units","degrees_north"))

          call check( nf90_def_var(ncid,'Mesh_face_lon',NF90_DOUBLE,dimid(2),varid))
            call check( nf90_put_att(ncid,varid,"standard_name","longitude"))
            call check( nf90_put_att(ncid,varid,"long_name","grid center longitude"))
            call check( nf90_put_att(ncid,varid,"units","degrees_east"))
            call check( nf90_put_att(ncid,varid,"bounds","Mesh_face_lonbnds"))

          call check( nf90_def_var(ncid,'Mesh_face_lat',NF90_DOUBLE,dimid(2),varid))
            call check( nf90_put_att(ncid,varid,"standard_name","latitude"))
            call check( nf90_put_att(ncid,varid,"long_name","grid center latitude"))
            call check( nf90_put_att(ncid,varid,"units","degrees_north"))
            call check( nf90_put_att(ncid,varid,"bounds","Mesh_face_latbnds"))

          call check( nf90_def_var(ncid,'Mesh_face_lonbnds',NF90_DOUBLE,(/dimid(3),dimid(2)/),varid))
          call check( nf90_def_var(ncid,'Mesh_face_latbnds',NF90_DOUBLE,(/dimid(3),dimid(2)/),varid))
        ENDIF

        IF (Parallel) THEN
          call check( nf90_def_var(ncid,'BulkElement_GlobalIndex',NF90_INT,dimid(2),varid))
            call check( nf90_put_att(ncid,varid,"mesh","Mesh"))
            call check( nf90_put_att(ncid,varid,"location","face"))

          call check( nf90_def_var(ncid,'Vertices_GlobalDOFs',NF90_INT,dimid(1),varid))
            call check( nf90_put_att(ncid,varid,"mesh","Mesh"))
            call check( nf90_put_att(ncid,varid,"location","node"))
        END IF

        call check( nf90_def_var(ncid,'BulkElement_Area',NF90_DOUBLE,dimid(2),varid))
          call check( nf90_put_att(ncid,varid,"mesh","Mesh"))
          call check( nf90_put_att(ncid,varid,"location","face"))

        call check( nf90_def_var(ncid,'time',NF90_DOUBLE,dimid(4),varid))

        ScalarFieldName = GetString( Params,'Scalar Field 1',ScalarsExist)
        Vari=1
        DO WHILE (ScalarsExist)
          Solution => VariableGet( Model % Mesh % Variables, TRIM(ScalarFieldName),ThisOnly=.TRUE.,UnFoundFatal=.TRUE.)
          VarType = Solution % TYPE
          SELECT CASE (VarType)
            CASE (Variable_on_nodes)
              call check( nf90_def_var(ncid,TRIM(Solution%Name),NF90_DOUBLE,(/dimid(1),dimid(4)/),varid))
                call check( nf90_put_att(ncid,varid,"mesh","Mesh"))
                call check( nf90_put_att(ncid,varid,"location","node"))
            CASE (Variable_on_elements)
              call check( nf90_def_var(ncid,TRIM(Solution%Name),NF90_DOUBLE,(/dimid(2),dimid(4)/),varid))
                call check( nf90_put_att(ncid,varid,"mesh","Mesh"))
                call check( nf90_put_att(ncid,varid,"location","face"))
            CASE DEFAULT
              CALL FATAL(Caller,"sorry don't know how to deal with this variable type")
          END SELECT

          Vari=Vari+1
          WRITE(Txt,'(A,I0)') 'Scalar Field ',Vari
          ScalarFieldName = GetString( Params,TRIM(Txt),ScalarsExist)
        END DO



        call check( nf90_enddef(ncid))
        call check( nf90_close(ncid))
      END  SUBROUTINE CreatNetcdfFile

      SUBROUTINE WriteMeshInfo(FName)
        IMPLICIT NONE
        CHARACTER(LEN=MAX_NAME_LEN),INTENT(IN) :: FName

        TYPE(Element_t),POINTER :: Element
        TYPE(Nodes_t),SAVE :: ElementNodes
        INTEGER, POINTER :: NodeIndexes(:)
        TYPE(GaussIntegrationPoints_t) :: IntegStuff
        REAL(KIND=dp) :: U,V,W,SqrtElementMetric
        REAL(KIND=dp),ALLOCATABLE :: Basis(:), dBasisdx(:,:)
        LOGICAL :: stat

        REAL(KIND=dp) :: xg,yg
        REAL(KIND=dp),DIMENSION(:),ALLOCATABLE :: x,y
        REAL(KIND=dp),DIMENSION(:),ALLOCATABLE :: NodeLon,FaceLon
        REAL(KIND=dp),DIMENSION(:),ALLOCATABLE :: NodeLat,FaceLat
        REAL(KIND=dp),DIMENSION(:,:),ALLOCATABLE :: LonBnds,LatBnds
        REAL(KIND=dp),DIMENSION(:),ALLOCATABLE :: area
        INTEGER,DIMENSION(:,:),ALLOCATABLE :: Indexes


        INTEGER :: i,ii,t,n,k
        INTEGER :: M
        INTEGER,DIMENSION(:),ALLOCATABLE :: Vertice
        INTEGER,DIMENSION(:),ALLOCATABLE :: GIndexes

        INTEGER :: VarId
        INTEGER :: ncid

        M = Model % MaxElementNodes
        ALLOCATE(Basis(M),dBasisdx(M,3))
        ALLOCATE(Vertice(NumberOfDofNodes),x(NumberOfDofNodes),y(NumberOfDofNodes))
        ALLOCATE(area(NumberOfElements),Indexes(MaxElementNodes,NumberOfElements),GIndexes(NumberOfElements))
        IF (ComputeLonLat) THEN
          ALLOCATE(NodeLon(NumberOfDofNodes),NodeLat(NumberOfDofNodes))
          ALLOCATE(FaceLon(NumberOfElements),FaceLat(NumberOfElements))
          ALLOCATE(LonBnds(MaxElementNodes,NumberOfElements),LatBnds(MaxElementNodes,NumberOfElements))
        END IF

        Indexes=Connect_Fill

        DO ii = 1, NumberOfDofNodes
          IF( NoPermutation ) THEN
            i = ii
          ELSE
            i = InvNodePerm(ii)
          END IF

          IF (Parallel) THEN
            Vertice(ii) = Mesh % ParallelInfo % GlobalDOFs(i)
          ELSE
            Vertice(ii) = i
          END IF
          x(ii) = Mesh%Nodes%x(i)
          y(ii) = Mesh%Nodes%y(i)

#ifdef HAVE_PROJ
          IF (ComputeLonLat) THEN
            coordp = pjuv_object(x(ii),y(ii))
            coordg = pj_inv(coordp, pj)
            NodeLon(ii) = coordg % u * pj_rad_to_deg
            NodeLat(ii) = coordg % v * pj_rad_to_deg
          ENDIF
#endif

        END DO

        area=0._dp
        t=0
        DO i = ElemFirst, ElemLast
          IF( .NOT. ActiveElem(i) ) CYCLE
          t=t+1

          Element => Model % Elements(i)
          IF (SaveLinear) THEN
            n = GetElementCorners( Element )
          ELSE
            n = GetElementNOFNodes(Element)
          END IF

          IF (Parallel) THEN
            GIndexes(t)=Element % GElementIndex
          ELSE
            GIndexes(t)=Element % ElementIndex
          ENDIF

          NodeIndexes => Element % NodeIndexes
          IF(IsClockwise( Mesh%Nodes%x(NodeIndexes(1:n)),Mesh%Nodes%y(NodeIndexes(1:n)),n)) THEN
            CALL FATAL(Caller,"Clock wise ordering ... implement reordering!")
          END IF

          IF( NoPermutation ) THEN
            Indexes(1:n,t) = NodeIndexes(1:n)
          ELSE
            Indexes(1:n,t) = NodePerm( NodeIndexes(1:n) )
          END IF

#ifdef HAVE_PROJ
          IF (ComputeLonLat) THEN
            ! element center
            xg=SUM(Mesh%Nodes%x(NodeIndexes(1:n)))/n
            yg=SUM(Mesh%Nodes%y(NodeIndexes(1:n)))/n
            coordp = pjuv_object(xg,yg)
            coordg = pj_inv(coordp, pj)
            FaceLon(t) =coordg % u * pj_rad_to_deg
            FaceLat(t) =coordg % v * pj_rad_to_deg

            DO k=1,n
              LonBnds(k,t)=NodeLon(Indexes(k,t))
              LatBnds(k,t)=NodeLat(Indexes(k,t))
            END DO
          ENDIF
#endif

          CALL GetElementNodes( ElementNodes, Element )

          IntegStuff = GaussPoints( Element )
          Do k=1,IntegStuff % n
            U = IntegStuff % u(k)
            V = IntegStuff % v(k)
            W = IntegStuff % w(k)
            stat = ElementInfo(Element,ElementNodes,U,V,W,SqrtElementMetric, &
                        Basis,dBasisdx )
            area(t)=area(t)+SqrtElementMetric*IntegStuff % s(k)
           End do

         END DO


        call check(nf90_open(trim(FName),NF90_WRITE,ncid))

        IF (Parallel) THEN
          call check(nf90_inq_varid(ncid,'Vertices_GlobalDOFs',VarId))
          call check(nf90_put_var(ncid,varid,Vertice(1:NumberOfDofNodes)))

          call check(nf90_inq_varid(ncid,'BulkElement_GlobalIndex',VarId))
          call check(nf90_put_var(ncid,varid,GIndexes))
        END IF

        call check(nf90_inq_varid(ncid,'Mesh_node_x',VarId))
        call check(nf90_put_var(ncid,varid,x(1:NumberOfDofNodes)))

        call check(nf90_inq_varid(ncid,'Mesh_node_y',VarId))
        call check(nf90_put_var(ncid,varid,y(1:NumberOfDofNodes)))

        call check(nf90_inq_varid(ncid,'Mesh_face_nodes',VarId))
        call check(nf90_put_var(ncid,VarId,Indexes))

        call check(nf90_inq_varid(ncid,'BulkElement_Area',VarId))
        call check(nf90_put_var(ncid,varid,area(1:NumberOfElements)))

        IF (ComputeLonLat) THEN
          call check(nf90_inq_varid(ncid,'Mesh_node_lon',VarId))
          call check(nf90_put_var(ncid,varid,NodeLon(1:NumberOfDofNodes)))

          call check(nf90_inq_varid(ncid,'Mesh_node_lat',VarId))
          call check(nf90_put_var(ncid,varid,NodeLat(1:NumberOfDofNodes)))

          call check(nf90_inq_varid(ncid,'Mesh_face_lat',VarId))
          call check(nf90_put_var(ncid,varid,FaceLat(1:NumberOfElements)))

          call check(nf90_inq_varid(ncid,'Mesh_face_lon',VarId))
          call check(nf90_put_var(ncid,varid,FaceLon(1:NumberOfElements)))

          call check(nf90_inq_varid(ncid,'Mesh_face_lonbnds',VarId))
          call check(nf90_put_var(ncid,varid,LonBnds))

          call check(nf90_inq_varid(ncid,'Mesh_face_latbnds',VarId))
          call check(nf90_put_var(ncid,varid,LatBnds))
        END IF

        call check(nf90_close(ncid))

        DEALLOCATE(Basis,dBasisdx)
        DEALLOCATE(Vertice,x,y)
        DEALLOCATE(area,Indexes,GIndexes)
        IF (ComputeLonLat) THEN
          DEALLOCATE(NodeLon,NodeLat)
          DEALLOCATE(FaceLon,FaceLat)
          DEALLOCATE(LonBnds,LatBnds)
        ENDIF

      END SUBROUTINE WriteMeshInfo

      SUBROUTINE WriteVariables(FName,nTime)
        IMPLICIT NONE
        CHARACTER(LEN=MAX_NAME_LEN),INTENT(IN) :: FName
        INTEGER,INTENT(IN) :: nTime

        TYPE(Element_t),POINTER :: Element
        TYPE(Variable_t),POINTER :: Solution
        INTEGER :: VarType
        REAL(KIND=dp) :: Time
        INTEGER :: ncid,VarId
        INTEGER :: ii,i,j,t,m

        REAL(KIND=dp),ALLOCATABLE :: NodeVar(:)
        REAL(KIND=dp),ALLOCATABLE :: EVar(:)

        INTEGER, POINTER :: Perm(:)
        REAL(KIND=dp),POINTER :: Values(:)

        CHARACTER(LEN=1024) :: Txt, ScalarFieldName
        INTEGER :: Vari
        LOGICAL ::ScalarsExist

        ALLOCATE(NodeVar(NumberOfDofNodes),EVar(NumberOfElements))

        call check(nf90_open(trim(FName),NF90_WRITE,ncid))
        
        Time = GetTime()
        call check(nf90_inq_varid(ncid,'time',VarId))
        call check(nf90_put_var(ncid,VarId,Time,start= (/nTime/)))

        ScalarFieldName = GetString( Params,'Scalar Field 1',ScalarsExist)
        Vari=1
        DO WHILE (ScalarsExist)
          Solution => VariableGet( Model % Mesh % Variables, TRIM(ScalarFieldName),ThisOnly=.TRUE.,UnFoundFatal=.TRUE.)
          VarType = Solution % TYPE
          Perm => Solution % Perm
          Values => Solution % Values
          SELECT CASE (VarType)
            CASE (Variable_on_nodes)
               DO ii=1,NumberOfDofNodes
                 IF( NoPermutation ) THEN
                   i = ii
                 ELSE
                   i = InvNodePerm(ii)
                 END IF
                 IF( ASSOCIATED( Perm ) ) THEN
                  j = Perm(i)
                 ELSE
                  j = i
                 END IF
                 IF (j==0) THEN
                  NodeVar(ii) = 0._dp !should make it fillValue?
                 ELSE
                  NodeVar(ii) = Values(j)
                 ENDIF
               END DO
               call check(nf90_inq_varid(ncid,TRIM(Solution%Name),VarId))
               call check(nf90_put_var(ncid,VarId,NodeVar,start= (/1,nTime/) ))
            CASE (Variable_on_elements)
               t=0
               DO i = ElemFirst, ElemLast
                 IF( .NOT. ActiveElem(i) ) CYCLE
                 t=t+1

                 Element => Model % Elements(i)
                 m = Element % ElementIndex
                 IF( ASSOCIATED( Perm ) ) THEN
                  IF( m>SIZE( Perm ) ) THEN
                    j = 0
                  ELSE
                    j = Perm(m)
                  ENDIF
                 ELSE
                  j = m
                 END IF
                 IF (j==0) THEN
                  EVar(t) = 0._dp !should make it fillValue?
                 ELSE
                  EVar(t) = Values(j)
                 ENDIF
               END DO
                 call check(nf90_inq_varid(ncid,TRIM(Solution%Name),VarId))
                 call check(nf90_put_var(ncid,VarId,EVar,start= (/1,nTime/) ))
            CASE DEFAULT
               CALL FATAL(Caller,"sorry don't know how to deal with this variable type")
          END SELECT

          Vari=Vari+1
          WRITE(Txt,'(A,I0)') 'Scalar Field ',Vari
          ScalarFieldName = GetString( Params,TRIM(Txt),ScalarsExist)
        END DO


        call check(nf90_close(ncid))
        
        DEALLOCATE(NodeVar,EVar)

      END SUBROUTINE WriteVariables


      subroutine check(status)
          integer, intent ( in) :: status
          CHARACTER(LEN=MAX_NAME_LEN) :: Message
          if(status /= nf90_noerr) then
            write(message,'(a)') trim(nf90_strerror(status))
            CALL FATAL(Caller,Message)
          end if
      end subroutine check

      FUNCTION IsClockwise(x,y,n) RESULT(clokwise)
      LOGICAL :: clokwise
      REAL(KIND=dp) :: x(n),y(n)
      INTEGER :: n
      INTEGER :: i,ind2
      REAL(KIND=dp) :: sarea

      sarea=0._dp
      Do i=1,n
        ind2=mod(i,n)+1
        sarea=sarea+(x(ind2)-x(i))*(y(ind2)+y(i))
      End do
      clokwise=(sarea.GT.0)

      END FUNCTION IsClockwise

      END SUBROUTINE NetcdfOutputSolver
