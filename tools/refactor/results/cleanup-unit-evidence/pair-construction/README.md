# B2/O1 candidate-source pair-construction probe

The five source files in source-manifest.json were evaluated from retained copies. The one native translation unit was freshly compiled at the displayed release flags; R facade functions used an earlier installed package only for unchanged dependencies. This is not an installed-package qualification.

61 independent checks in three tests passed with zero failures, warnings, errors or skips. They cover inclusive unsorted and wide-coordinate selection, safe public midpoint addition, public row ordering, trans Cartesian-product ordering, both caller table classes and serialization/nonmutation. Six separately executed old/new trans cases compare entire returned objects including base expand.grid out.attrs using independent input copies.

The first run passed all predicate/cis checks but seven trans comparisons expected a plain independently built data.frame without existing expand.grid out.attrs. expected-attributes-attempt.log preserves that test-only mismatch. The corrected tests compare independent columns and class; compare-trans.R separately checks complete original/candidate attributes. No production behavior was changed in response.

The original signed-overflow/missing-coordinate native cases were not executed. Source preconditions protect unrepresentable native coordinates without inventing a biological missing-coordinate policy. Public formal arguments are unchanged. Final combined immutable installed, native and consumer qualification remains required.
