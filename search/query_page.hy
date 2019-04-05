(import [pyhtml [*]])

(import [search.template_common [paths-to-static]])

(require [hy.contrib.walk [let]])

(require [search.template_common [build-page]])

(defn form-required []
  """Build the portion of the query form that must be filled
out for a search.
"""
  (fieldset
    (legend "Search")

    (let [field-name "genome"]
      ((div :class_ "pure-control-group")
       ((label :for_ field-name)
        "Genome: ")

       ((select :id field-name :name field-name)
        ((option :value "dm6")
         "D. melanogaster"))))

    (let [field-name "sequences"]
      ((div :class_ "pure-control-group")
       ((label :for_ field-name)
        "Sequence(s): "
        ((span :class_ "pure-form-message")
         "One sequence per line"))

       (textarea
         :id field-name
         :name field-name
         :cols "30"
         :rows "2")))

    (let [field-name "submit"]
      ((div :class_ "pure-controls")
       ((button :id field-name :type "submit" :class_ "pure-button pure-button-primary")
        "Submit")))))

(defn build-query-page []
  "Build and return the query page's Django template."
  (let [generated
        (build-page
          "Query"
          ((form :class_ "pure-form pure-form-aligned")
           (form-required)))]
    (paths-to-static generated)))
