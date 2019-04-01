// Query {
//   genome: String,
//   queries: [String, ….], //sequences
// }

function buildQuery() {
  let elements = {
    "genome": document.getElementById("genome"),
    "sequences": document.getElementById("sequences")
  };

  let res = {
    "genome": null,
    "sequences": null
  }

  for (key in elements) {
    if (elements[key] == null) {
      // TODO: Properly throw an exception.
      alert(`elements[${key}] is null`);
    }
  }

  res["genome"] = elements["genome"].innerText;
  res["sequences"] = elements["sequences"].value
    .split("\n")
    .filter((sequence) => { return sequence.length > 0; });

  return res;
}

document.getElementById("submit").addEventListener("click", (evt) => {
  evt.preventDefault();
  console.log(buildQuery());
});
