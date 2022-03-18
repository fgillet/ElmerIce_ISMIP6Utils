      FUNCTION IcyMask(Model,nodenumber,VarIn) RESULT(VarOut)
      USE DefUtils
      implicit none
      !-----------------
      TYPE(Model_t) :: Model
      INTEGER :: nodenumber
      REAL(kind=dp) :: VarIn,VarOut

      TYPE(Element_t), POINTER :: Element
      TYPE(Variable_t),POINTER :: Var
      TYPE(ValueList_t), POINTER :: Material

      INTEGER :: n

      REAL(KIND=dp),ALLOCATABLE,SAVE :: NodalVar(:),MinH(:)
      LOGICAL,SAVE :: FirstTime=.TRUE.

      IF (FirstTime) THEN
         n= Model % MaxElementNodes
         allocate(NodalVar(n),MinH(n))
         FirstTime=.FALSE.
      ENDIF

      !Get current element
      Element => GetCurrentElement()
      n = GetElementNOFNodes()
      Material => GetMaterial(Element)

      !Get groundedMask
      Var =>  VariableGet(Model % Variables, 'groundedMask',UnFoundFatal=.TRUE.)

      NodalVar(1:n)=Var % values ( Var % Perm (Element % NodeIndexes(1:n)))
      IF (ALL(NodalVar(1:n).GT.-0.5)) THEN
       VarOut = +1
      ELSE IF (ALL(NodalVar(1:n).LT.-0.5)) THEN
       VarOut = +2
      ELSE
       VarOut = 0
      ENDIF

      !Get Thickness
      Var =>  VariableGet(Model % Variables, 'h',UnFoundFatal=.TRUE.)
      NodalVar(1:n)=Var % values ( Var % Perm (Element % NodeIndexes(1:n)))

      MinH=0._dp
      MinH(1:n) =  ListGetReal( Material,'Min H',n,Element % NodeIndexes,UnFoundFatal=.TRUE.)

      IF (ALL((NodalVar(1:n)-MinH(1:n)).LE.0._dp)) THEN
       VarOut = -1
      END IF

      End FUNCTION IcyMask

