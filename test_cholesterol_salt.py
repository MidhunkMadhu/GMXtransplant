import unittest

from cholesterol import compare_salt_concentration


class CholesterolSaltComparisonTests(unittest.TestCase):
    def compare(self, final_pairs):
        return compare_salt_concentration(
            {"SOD": 143, "CLA": 150},
            {"SOD": final_pairs, "CLA": final_pairs + 7},
            "SOD",
            "CLA",
            0.15,
            0.05,
        )

    def test_unchanged_reference_concentration(self):
        result = self.compare(143)
        self.assertEqual(result["reference_salt_pairs"], 143)
        self.assertEqual(result["final_salt_pairs"], 143)
        self.assertAlmostEqual(result["estimated_final_concentration_molar"], 0.15)
        self.assertTrue(result["concentration_within_tolerance"])

    def test_seven_pair_loss_is_within_five_percent(self):
        result = self.compare(136)
        self.assertLess(result["concentration_deviation_fraction"], 0.05)
        self.assertTrue(result["concentration_within_tolerance"])

    def test_eight_pair_loss_exceeds_five_percent(self):
        result = self.compare(135)
        self.assertGreater(result["concentration_deviation_fraction"], 0.05)
        self.assertFalse(result["concentration_within_tolerance"])

    def test_missing_reference_pairs_is_reportable(self):
        result = compare_salt_concentration(
            {"SOD": 0, "CLA": 7},
            {"SOD": 0, "CLA": 7},
            "SOD",
            "CLA",
            0.15,
            0.05,
        )
        self.assertIsNone(result["estimated_final_concentration_molar"])
        self.assertIsNone(result["concentration_within_tolerance"])


if __name__ == "__main__":
    unittest.main()
