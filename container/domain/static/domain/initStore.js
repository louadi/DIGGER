"use strict"
;(() => {
  const PREFIX = "nease/"

  const storage = {
    set: (key, value) => {
      try {
        if (value == null) {
          localStorage.removeItem(key)
        } else {
          localStorage.setItem(key, JSON.stringify(value))
        }
      } catch {}
    },
    get: key => {
      const value = localStorage.getItem(key)
      if (value == null) return null
      try {
        return JSON.parse(value)
      } catch {
        localStorage.removeItem(key)
        return null
      }
    },
  }

  const getExpiration = expiresInDays => {
    if (expiresInDays == null) return null
    const expiresInMs = expiresInDays * 24 * 60 * 60 * 1000
    return Date.now() + expiresInMs
  }

  const validateOptions = options => {
    if (options == null || typeof options !== "object") {
      throw new Error("[initStore]: options must be passed")
    }
    const { key, expiresInDays, name } = options
    if (key == null || typeof key !== "string") {
      throw new Error(
        "[initStore]: options.key is required and must be a string"
      )
    }
    if (expiresInDays && typeof expiresInDays !== "number") {
      throw new Error("[initStore]: options.expiresInDays must be a number")
    }
    if (name && typeof name !== "string") {
      throw new Error("[initStore]: options.name must be a string")
    }
    Object.keys(options).forEach(key => {
      if (!["key", "expiresInDays", "name"].includes(key)) {
        throw new Error(`[initStore]: unknown option ${key}`)
      }
    })

    return options
  }

  /** Function to create a store to persist a value in the localStorage
   *
   * @param {Object} options
   * @param {string} options.key
   * @param {number | undefined} options.expiresInDays
   * @param {string | undefined} options.name
   *
   * @returns {Object} store with functions to set and get the value
   **/
  const initStore = options => {
    console.log("Initializing store with options", options)
    const { key, expiresInDays, name } = validateOptions(options)
    const prefixedKey = PREFIX + key
    const prefixedName = name ? PREFIX + name : null

    let timeout = null
    const setExpiration = expiresAt => {
      clearTimeout(timeout)
      if (expiresAt != null) {
        console.log("Setting timeout for expiration", expiresAt)
        timeout = setTimeout(() => {
          console.log("Removing item from store")
          storage.set(prefixedKey, null)
          storage.set(prefixedName, null)
        }, expiresAt - Date.now())
      }
    }

    const { expiresAt } = storage.get(prefixedKey) ?? {}

    if (expiresAt != null) {
      if (expiresAt < Date.now()) {
        storage.set(prefixedKey, null)
        storage.set(prefixedName, null)
      } else {
        setExpiration(expiresAt)
      }
    }

    return {
      get: () => storage.get(prefixedKey)?.value ?? null,
      set: value => {
        const newValue =
          value == null
            ? null
            : { value, expiresAt: getExpiration(expiresInDays), name }
        setExpiration(newValue?.expiresAt)
        storage.set(prefixedKey, newValue)
      },
    }
  }
  window.initStore = initStore
})()

/** Usage Examples: **/

/* Store that persists a user id

const userStore = window.initStore({ key: "user-id" })
userStore.set("177013")
const userId = userStore.get()
 
*/

/* Store that expires after 5 days

const userStore = window.initStore({ key: "user-id", expiresInDays: 5 })
userStore.set("177013")
// value will be stored for 5 days or until setting the value again
 
*/
