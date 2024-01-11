from utilities.parse import xyz_parse
from Excpetions import TAFFIException

class StructureException(TAFFIException):
    pass

class StructureBase:
    """
    base class for almost everything, holds elements and geometry
    """
    def __init__(self):
        self._elements = None
        self._geometry = None
        self._q_tot = 0
        # often times the structure comes from a xyz,
        # keeping the xyz name is useful for debugging
        self.xyz = None

    # base class needs to be protected
    @property
    def elements(self):
        if self._elements is not None:
            return self._elements
        else:
            raise StructureException("Elements not parsed")
    @elements.setter
    def elements(self, value):
        self._elements = value

    @property
    def geometry(self):
        if self._geometry is not None:
            return self._geometry
        else:
            raise StructureException("Geometry not parsed")

    @geometry.setter
    def geometry(self, value):
        self._geometry = value

    @property
    def q_tot(self):
        return self._q_tot

    @q_tot.setter
    def q_tot(self, value):
        self._q_tot = value

    def parse_from_xyz(self, xyz, q_opt=False):
        self.xyz = xyz
        if q_opt:
            self.elements, self.geometry, self.q_tot = xyz_parse(xyz, q_opt=q_opt)
        else:
            self.elements, self.geometry = xyz_parse(xyz, q_opt=q_opt)

    def parse_data(self, **kwargs):
        elements = kwargs.get('elements', None)
        geometry = kwargs.get('geometry', None)
        q_tot = kwargs.get('q_tot', None)
        if elements is not None:
            self.elements = elements
        if geometry is not None:
            self.geometry = geometry
        if q_tot is not None:
            self.q_tot = q_tot

def main():
    struc = StructureBase()
    struc.parse_from_xyz('test.xyz')
    test_dict = {'option2':1, 'option3':2}


if __name__ == '__main__':
    main()
