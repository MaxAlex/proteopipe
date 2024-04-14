import bisect

class SortedList:
    def __init__(self, key_func, initial_items=None):
        """Initialize with a key function to sort & search by and an optional list of initial items."""
        self.key_func = key_func
        if initial_items is None:
            self._keys = []
            self._items = []
        else:
            # Sort the initial items using the key function
            sorted_items = sorted(initial_items, key=self.key_func)
            self._keys = [self.key_func(item) for item in sorted_items]
            self._items = sorted_items

    def __getitem__(self, key):
        """Find an item by its key using binary search.
        
        Return the closest item if found, otherwise raise an error if the list is empty.
        """
        if not self._items:
            raise ValueError("The list is empty.")
        
        index = bisect.bisect_left(self._keys, key)
        # if the exact key is not found, we provide the closest item
        if index == len(self._keys) or (index > 0 and key - self._keys[index-1] < self._keys[index] - key):
            index -= 1
        return self._items[index]

    def __str__(self):
        return str(self._items)


mod_by_mass = {57.0214645:'Carbamidomethyl',
               15.9949155:'Oxidation',
               42.0105655:'Acetyl',
               79.9663005:'Phospho',
               229.16295: 'TMT6plex'}
mod_masses = SortedList(lambda x: x[0], mod_by_mass.items())


def identify_mod_by_mass(modmass, da_tol = 0.01):
    # Get nearest mod mass, to a tolerance of +- 0.01 Da,
    # and return corresponding modification name

    near_mass, near_mod = mod_masses[modmass]
    if abs(near_mass - modmass) < da_tol:
        return near_mod
    else:
        raise Exception("Unknown modification mass %s" % modmass)


if __name__ == '__main__':
    test_cases = [
        (57.0214645, 'Carbamidomethyl'),
        (57.0214647, 'Carbamidomethyl'),
        (15.9949160, 'Oxidation'),
        (42.0105600, 'Acetyl'),
        (80.9663000, 'Phospho'),
        (229.16300, 'TMT6plex'),
    ]
    
    for mass, expected_mod in test_cases:
        matched_mod = mod_masses[mass][1]
        assert matched_mod == expected_mod, f"Expected: {expected_mod}, Got: {matched_mod}"
    print("Tests passed")
